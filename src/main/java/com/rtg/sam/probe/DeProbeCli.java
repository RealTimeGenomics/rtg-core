/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
package com.rtg.sam.probe;

import static com.rtg.launcher.CommonFlags.MIN_READ_LENGTH;
import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.io.UncheckedIOException;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;

import com.rtg.bed.BedReader;
import com.rtg.bed.BedRecord;
import com.rtg.bed.BedWriter;
import com.rtg.launcher.CommandLineFiles;
import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.LoggedCli;
import com.rtg.sam.ExtraSoftClip;
import com.rtg.sam.MinLengthFilter;
import com.rtg.sam.RecordIterator;
import com.rtg.sam.SamFilterParams;
import com.rtg.sam.SamOutput;
import com.rtg.sam.SamReadingContext;
import com.rtg.sam.SamRecordPopulator;
import com.rtg.sam.SamUtils;
import com.rtg.sam.SmartSamWriter;
import com.rtg.sam.ThreadedMultifileIterator;
import com.rtg.util.Counter;
import com.rtg.util.MathUtils;
import com.rtg.util.SingletonPopulatorFactory;
import com.rtg.util.StringUtils;
import com.rtg.util.TextTable;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Flag;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.intervals.RangeList;
import com.rtg.util.intervals.SimpleRangeMeta;
import com.rtg.util.intervals.RangeMeta;
import com.rtg.util.intervals.ReferenceRanges;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.LogStream;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;

/**
 * Remove probes from mapped reads using a bed file to describe probes.
 * <ul>
 *   <li>Assumes probe is on first read in sequencing</li>
 *   <li>Assumes 6th column is +/- representing whether probe is on forward or reverse strand</li>
 *   <li>Only updates records that are mapped, and marked as first in sequencing</li>
 *   <li>Update is made by removing values as appropriate from: the read string, the quality string, and the cigar. As well as updating the alignment start position if required.</li>
 *   <li>Mapped first arm reads that do not match a known probe will receive an annotation (XS:Z:failed)</li>
 * </ul>
 */
public class DeProbeCli extends LoggedCli {

  static final String ATTRIBUTE_STRIPPED_PROBE = "XP";
  static final String ATTRIBUTE_STRIP_STATUS = "XS";

  private static final String PROBE_DROPPING = "probe-specific-filter";
  private static final String PROBE_BED = "probe-bed";
  private static final String TOLERANCE_FLAG = "tolerance";
  private static final String EXTRA_SOFT_CLIP_FLAG = "extra-soft-clip";

  static final String ALIGNMENT_FILE_NAME = "alignments_deprobed.bam";
  static final String PROBE_OFFSET_TABLE_FILE = "probe_offsets.tsv";
  static final String CIGAR_OP_TABLE_FILE = "cigar_ops.tsv";
  static final String PROBE_SUMMARY_FILE = "probe_summary.tsv";
  static final String ON_TARGET_SUMMARY_FILE = "on_target.tsv";
  static final String COUNTS_SUFFIX = "_strand_probe_counts.bed.gz";
  static final String POS_COUNTS_NAME = "positive" + COUNTS_SUFFIX;
  static final String NEG_COUNTS_NAME = "negative" + COUNTS_SUFFIX;


  private ReferenceRanges<ProbeCounter> mPosRanges;
  private ReferenceRanges<ProbeCounter> mNegRanges;
  private long mTotalMappedReads;
  private long mTotalMappedPos;
  private long mTotalMappedNeg;
  private long mTotalStrippedPos;
  private long mTotalStrippedNeg;
  private long mTotalStrippedReads;


  @Override
  public String moduleName() {
    return "deprobe";
  }

  @Override
  public String description() {
    return "trim strand-specific probes from mapped SAM/BAM alignments";
  }

  @Override
  protected void initFlags() {
    mFlags.setDescription(StringUtils.sentencify(description()));
    CommonFlagCategories.setCategories(mFlags);
    final Flag<File> inFlag = mFlags.registerRequired(File.class, CommonFlags.FILE, "SAM/BAM format files containing coordinate-sorted alignments")
      .setCategory(INPUT_OUTPUT).setMinCount(0).setMaxCount(Integer.MAX_VALUE);
    final Flag<File> listFlag = mFlags.registerOptional('I', CommonFlags.INPUT_LIST_FLAG, File.class, CommonFlags.FILE, "file containing a list of SAM/BAM format files (1 per line) containing coordinate-sorted alignments").setCategory(INPUT_OUTPUT);
    CommonFlags.initOutputDirFlag(mFlags);
    mFlags.registerRequired('b', PROBE_BED, File.class, CommonFlags.FILE, "BED file specifying each probe location and strand").setCategory(INPUT_OUTPUT);
    mFlags.registerOptional(TOLERANCE_FLAG, Integer.class, CommonFlags.INT, "start position tolerance for probe matching", 5).setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    mFlags.registerOptional(EXTRA_SOFT_CLIP_FLAG, "if set, add extra soft-clipping where mismatches occur at the end of reads").setCategory(CommonFlagCategories.FILTERING);
    mFlags.registerOptional(PROBE_DROPPING, "if set, apply probe-specific minimum read length filter as supplied in probe bed column 5").setCategory(CommonFlagCategories.FILTERING);
    CommonFlags.initMinReadLength(mFlags);
    mFlags.addRequiredSet(inFlag);
    mFlags.addRequiredSet(listFlag);

    mFlags.setValidator(flags -> CommonFlags.validateOutputDirectory(flags)
      && CommonFlags.checkFileList(flags, CommonFlags.INPUT_LIST_FLAG, null, Integer.MAX_VALUE)
      && CommonFlags.validateInputFile(flags, PROBE_BED)
      && flags.checkInRange(TOLERANCE_FLAG, 0, Integer.MAX_VALUE));
  }

  @Override
  protected File outputDirectory() {
    return (File) mFlags.getValue(CommonFlags.OUTPUT_FLAG);
  }

  @Override
  protected int mainExec(OutputStream out, LogStream log) throws IOException {
    final int tolerance = (Integer) mFlags.getValue(TOLERANCE_FLAG);
    final boolean extraSoftClip = mFlags.isSet(EXTRA_SOFT_CLIP_FLAG);
    final Collection<File> inputFiles = new CommandLineFiles(CommonFlags.INPUT_LIST_FLAG, null, CommandLineFiles.EXISTS, CommandLineFiles.NOT_DIRECTORY).getFileList(mFlags);
    final PosChecker posChecker = new PosChecker(tolerance);
    final NegChecker negChecker = new NegChecker(tolerance);
    final File bedFile = (File) mFlags.getValue(PROBE_BED);
    final int threads = CommonFlags.parseIOThreads(null); // Just use the default
    final int minBases = (int) mFlags.getValue(MIN_READ_LENGTH);
    final boolean probeBasedDropping = mFlags.isSet(PROBE_DROPPING);

    mTotalMappedReads = 0;
    mTotalMappedPos = 0;
    mTotalMappedNeg = 0;
    mTotalStrippedPos = 0;
    mTotalStrippedNeg = 0;
    mTotalStrippedReads = 0;

    int totalRecords12 = 0;
    int softClippedRecords12 = 0;
    int shortReads = 0;

    boolean warnedNoReadGroup = false;
    long totalUnmappedReads = 0;
    long totalRecords = 0;
    long totalAmbiguousRecords = 0;
    final HashMap<String, ReadStatus> ambiguousReads = new HashMap<>();
    String currentSequence = null;
    final SamReadingContext c = new SamReadingContext(inputFiles, threads, new SamFilterParams.SamFilterParamsBuilder().create(), SamUtils.getUberHeader(inputFiles), null);
    try (RecordIterator<SAMRecord> reader = new ThreadedMultifileIterator<>(c, new SingletonPopulatorFactory<>(new SamRecordPopulator()))) {
      loadProbeBed(bedFile);
      final File outputFile = new File(outputDirectory(), ALIGNMENT_FILE_NAME);
      try (SamOutput samOutput = SamOutput.getSamOutput(outputFile, out, reader.header(), !mFlags.isSet(CommonFlags.NO_GZIP), true, null)) {
        try (SmartSamWriter writer = new SmartSamWriter(samOutput.getWriter())) {
          while (reader.hasNext()) {
            final SAMRecord record = reader.next();
            if (!record.getReferenceName().equals(currentSequence)) {
              Diagnostic.progress("Processing: " + record.getReferenceName());
              currentSequence = record.getReferenceName();
            }
            // Update read name to include rgid to handle reads mapped as IDs rather than full names
            final SAMReadGroupRecord readGroup = record.getReadGroup();
            if (readGroup != null && readGroup.getId() != null) {
              record.setReadName(readGroup.getId() + "-" + record.getReadName());
            } else if (!warnedNoReadGroup) {
              Diagnostic.warning("Encountered alignment without read group information. Statistics may not be correct if mappings to not include full read names");
              warnedNoReadGroup = true;
            }
            if (record.getFirstOfPairFlag()) {
              ++totalRecords;
              final int ih = MathUtils.unboxNatural(SamUtils.getNHOrIH(record));
              final boolean unique = 1 == ih || !record.isSecondaryAlignment();
              final boolean negative = record.getReadNegativeStrandFlag();
              final boolean mapped = !record.getReadUnmappedFlag();
              boolean stripped = false;
              if (mapped) {
                ProbeCounter hitProbe = checkList(mPosRanges.get(record.getReferenceName()), record, null, tolerance, posChecker);
                if (hitProbe == null) {
                  hitProbe = checkList(mNegRanges.get(record.getReferenceName()), record, null, tolerance, negChecker);
                }
                stripped = hitProbe != null;
                if (!stripped) {
                  record.setAttribute(ATTRIBUTE_STRIP_STATUS, "failed"); // Off-target alignment
                } else if (probeBasedDropping && record.getReadBases().length < hitProbe.getMinReadLength()) {
                  SamUtils.convertToUnmapped(record);
                  shortReads++;
                } else {
                  // Sometimes after stripping we're just left with soft-clipped bases. Convert these to unmapped
                  if (!record.getCigar().getCigarElements().stream().anyMatch(e -> e.getOperator().isIndelOrSkippedRegion() || e.getOperator().isAlignment())) {
                    //Diagnostic.warning("Encountered trimmed alignment without a real operator: " + record.getAlignmentStart() + "-" + record.getAlignmentEnd() + " " + record.getSAMString());
                    SamUtils.convertToUnmapped(record);
                    shortReads++;
                  }
                }
              }

              if (!mapped) {
                ++totalUnmappedReads;
              } else if (unique) {
                addMapped(stripped, negative);
              } else {
                // If not uniquely mapped, accumulate and delay addition to the stats until the end
                ++totalAmbiguousRecords;
                final String name = record.getReadName();
                ReadStatus s = ambiguousReads.get(name);
                if (s == null) {
                  s = new ReadStatus(name, stripped, negative);
                  ambiguousReads.put(name, s);
                } else if (stripped && !s.stripped()) {
                  s.setStatus(stripped, negative);
                }
              }
            }
            totalRecords12++;
            if (minBases > 0 && MinLengthFilter.filterShortReads(record, minBases)) {
              shortReads++;
            }

            if (extraSoftClip && ExtraSoftClip.addSoftClip(record)) {
              softClippedRecords12++;
            }

            if (totalRecords12 % 1000000 == 0) {
              logStats(extraSoftClip, minBases > 0 || probeBasedDropping, totalRecords12, softClippedRecords12, shortReads);
            }

            writer.addRecord(record);
          }

          logStats(extraSoftClip, minBases > 0 || probeBasedDropping, totalRecords12, softClippedRecords12, shortReads);
          if (writer.getDuplicateCount() > 0) {
            Diagnostic.warning(writer.getDuplicateCount() + " duplicate records were dropped during output.");
          }
          Diagnostic.userLog("Reordering buffer used capacity of " + writer.getMaxCapacityUsed() + " records");
          Diagnostic.progress("Closing SAM writer");
        }
      }
    }

    final long totalReads = mTotalMappedReads + totalUnmappedReads + ambiguousReads.size();
    Diagnostic.userLog("Total R1+R2 records: " + totalRecords12);
    Diagnostic.userLog("Total R1+R2 records soft-clipped: " + softClippedRecords12);
    Diagnostic.userLog("Total R1 records: " + totalRecords);
    Diagnostic.userLog("Total R1 reads: " + totalReads);
    Diagnostic.userLog("Total R1 unmapped reads (records): " + totalUnmappedReads);
    Diagnostic.userLog("Total R1 uniquely mapped reads (records): " + mTotalMappedReads);
    Diagnostic.userLog("Total R1 ambiguously mapped records: " + totalAmbiguousRecords);
    Diagnostic.userLog("Total R1 ambiguously mapped reads: " + ambiguousReads.size());

    // Now merge in counts of any ambiguously mapped reads
    for (ReadStatus s : ambiguousReads.values()) {
      addMapped(s.stripped(), s.negative());
    }

    reportOffsetSummary(posChecker, negChecker, tolerance);
    reportCigarSummary(posChecker, negChecker);
    reportProbeSummary(posChecker, negChecker);
    reportOnTargetSummary(out, totalReads);
    writeStrandFile(new File(outputDirectory(), POS_COUNTS_NAME), mPosRanges);
    writeStrandFile(new File(outputDirectory(), NEG_COUNTS_NAME), mNegRanges);
    return 0;
  }

  private void logStats(boolean extraSoftClip, boolean minBases, int totalRecords12, int softClippedRecords12, int shortReads) {
    final StringBuilder logString = new StringBuilder();
    if (minBases) {
      logString.append(String.format(", Short reads: %d (%.2f%%)", shortReads, 100.0 * shortReads / totalRecords12));
    }
    if (extraSoftClip) {
      logString.append(String.format(", Soft-clipped reads: %d (%.2f%%)", softClippedRecords12, 100.0 * softClippedRecords12 / totalRecords12));
    }
    Diagnostic.developerLog(String.format("Records: %d", totalRecords12) + logString.toString());
  }

  private void writeStrandFile(File file, ReferenceRanges<ProbeCounter> posRanges) throws IOException {
    try (OutputStream posFile = FileUtils.createOutputStream(file)) {
      writeProbeCounts(posFile, posRanges);
    }
  }

  private void reportOnTargetSummary(OutputStream out, long totalReads) throws IOException {
    Diagnostic.userLog("ON TARGET RATES");
    final TextTable onTargetSummary = new TextTable();
    onTargetSummary.addHeaderRow("group", "total", "mapped", "on target", "fraction of total", "fraction of mapped");
    onTargetSummary.addSeparator();
    onTargetSummary.addRow("Reads", Long.toString(totalReads), Long.toString(mTotalMappedReads), Long.toString(mTotalStrippedReads),
      String.format("%.4f", (double) mTotalStrippedReads / (double) totalReads),
      String.format("%.4f", (double) mTotalStrippedReads / (double) mTotalMappedReads));
    Diagnostic.userLog(onTargetSummary.toString());
    FileUtils.stringToFile(onTargetSummary.getAsTsv(), new File(outputDirectory(), ON_TARGET_SUMMARY_FILE));

    try (PrintStream summaryOut = new PrintStream(FileUtils.createTeedOutputStream(FileUtils.createOutputStream(new File(outputDirectory(), CommonFlags.SUMMARY_FILE)), out))) {
      summaryOut.println(onTargetSummary);
    }
  }

  private void reportOffsetSummary(PosChecker posChecker, NegChecker negChecker, int tolerance) throws IOException {
    Diagnostic.userLog("PROBE OFFSETS");
    final TextTable offsetSummary = new TextTable();
    offsetSummary.addHeaderRow("delta", "+", "-");
    offsetSummary.addSeparator();
    for (int i = 0; i < posChecker.mPosDiffStats.length; ++i) {
      offsetSummary.addRow(Integer.toString(i - tolerance), Integer.toString(posChecker.mPosDiffStats[i]), Integer.toString(negChecker.mPosDiffStats[tolerance * 2 - i]));
    }
    Diagnostic.userLog(offsetSummary.toString());
    FileUtils.stringToFile(offsetSummary.getAsTsv(), new File(outputDirectory(), PROBE_OFFSET_TABLE_FILE));
  }

  private void reportCigarSummary(PosChecker posChecker, NegChecker negChecker) throws IOException {
    Diagnostic.userLog("CIGAR OPERATIONS WITHIN PROBE");
    final TextTable cigarSummary = new TextTable();
    cigarSummary.addHeaderRow("strand", "length", "X", "I", "D", "S");
    cigarSummary.addSeparator();
    for (int i = 0; i < PositionAndStrandChecker.MAX_OP_LEN; ++i) {
      addCigarRow(cigarSummary, "+", posChecker, i);
      addCigarRow(cigarSummary, "-", negChecker, i);
    }
    Diagnostic.userLog(cigarSummary.toString());
    FileUtils.stringToFile(cigarSummary.getAsTsv(), new File(outputDirectory(), CIGAR_OP_TABLE_FILE));
  }

  private void reportProbeSummary(PosChecker posChecker, NegChecker negChecker) throws IOException {
    final long totalstripped = mTotalStrippedPos + mTotalStrippedNeg;
    final long total = mTotalMappedPos + mTotalMappedNeg;
    final TextTable summary = new TextTable();
    summary.addHeaderRow("strand", "alignments", "stripped", "fraction", "nt/read");
    summary.addSeparator();
    summary.addRow("+", Long.toString(mTotalMappedPos), Long.toString(mTotalStrippedPos),
      String.format("%.4f", (double) mTotalStrippedPos / mTotalMappedPos),
      String.format("%.1f", (double) posChecker.mBasesTrimmed / mTotalStrippedPos));
    summary.addRow("-", Long.toString(mTotalMappedNeg), Long.toString(mTotalStrippedNeg),
      String.format("%.4f", (double) mTotalStrippedNeg / mTotalMappedNeg),
      String.format("%.1f", (double) negChecker.mBasesTrimmed / mTotalStrippedNeg));
    summary.addRow("Both", Long.toString(total), Long.toString(totalstripped),
      String.format("%.4f", (double) totalstripped / total),
      String.format("%.1f", (double) (posChecker.mBasesTrimmed + negChecker.mBasesTrimmed) / totalstripped));
    FileUtils.stringToFile(summary.getAsTsv(), new File(outputDirectory(), PROBE_SUMMARY_FILE));
  }
  private void writeProbeCounts(OutputStream output, ReferenceRanges<ProbeCounter> posRanges) throws IOException {
    try (BedWriter bedWriter = new BedWriter(output)) {
      bedWriter.writeComment("chrom\tstart\tend\tname\tscore\tstrand\tcount");
      for (String sequenceName : posRanges.sequenceNames()) {
        final RangeList<ProbeCounter> rangeList = posRanges.get(sequenceName);
        rangeList.getRangeList().stream().map(RangeList.RangeView::getEnclosingRanges).flatMap(Collection::stream).forEach(
          data -> {
            try {
              bedWriter.write(toBedRecord(sequenceName, data));
            } catch (IOException e) {
              throw new UncheckedIOException(e);
            }
          }
        );
      }
    } catch (UncheckedIOException e) {
      throw e.getCause();
    }
  }

  BedRecord toBedRecord(String sequenceName, RangeMeta<ProbeCounter> data) {
    final ProbeCounter counter = data.getMeta();
    final List<String> annotations = counter.getAnnotations();
    final String[] annotationArray = new String[annotations.size() + 1];
    int i = 0;
    for (String annotation : annotations) {
      annotationArray[i++] = annotation;
    }
    annotationArray[i] = Integer.toString(counter.count());
    return new BedRecord(sequenceName, data.getStart(), data.getEnd(), annotationArray);
  }

  private void addMapped(boolean stripped, boolean negative) {
    ++mTotalMappedReads;
    if (stripped) {
      ++mTotalStrippedReads;
    }
    if (negative) {
      ++mTotalMappedNeg;
      if (stripped) {
        ++mTotalStrippedNeg;
      }
    } else {
      ++mTotalMappedPos;
      if (stripped) {
        ++mTotalStrippedPos;
      }
    }
  }

  private static class ProbeCounter extends Counter {
    private final List<String> mAnnotations;
    private final int mMinReadLength;

    ProbeCounter(List<String> annotations) {
      this.mAnnotations = annotations;
      int minReadLength = -1;
      if (annotations.size() > 1) {
        try {
          minReadLength = Integer.parseInt(annotations.get(1));
        } catch (NumberFormatException nfe) {
          throw new NoTalkbackSlimException("Expected Integer minimum read length in column 5 for probe " + annotations.get(0) + " was " + annotations.get(1));
        }
      }
      mMinReadLength = minReadLength;

    }

    public List<String> getAnnotations() {
      return mAnnotations;
    }

    public int getMinReadLength() {
      return mMinReadLength;
    }
  }


  private void addCigarRow(TextTable ctype, String strand, PositionAndStrandChecker checker, int i) {
    final String len = Integer.toString(i + 1) + (i == PositionAndStrandChecker.MAX_OP_LEN - 1 ? "+" : "");
    ctype.addRow(strand, len, Integer.toString(checker.mMismatchStats[i]), Integer.toString(checker.mInsertStats[i]), Integer.toString(checker.mDeletionStats[i]), Integer.toString(checker.mSoftClipStats[i]));
  }

  private void loadProbeBed(File bedFile) throws IOException {
    final ReferenceRanges.Accumulator<ProbeCounter> posRangesAccum = new ReferenceRanges.Accumulator<>();
    final ReferenceRanges.Accumulator<ProbeCounter> negRangesAccum = new ReferenceRanges.Accumulator<>();
    try (BedReader bedReader = BedReader.openBedReader(null, bedFile, 3)) {
      while (bedReader.hasNext()) {
        final BedRecord record = bedReader.next();
        final RangeMeta<ProbeCounter> rangeData = new SimpleRangeMeta<>(record.getStart(), record.getEnd(), new ProbeCounter(Arrays.asList(record.getAnnotations())));
        if ("+".equals(record.getAnnotations()[2])) {
          posRangesAccum.addRangeData(record.getSequenceName(), rangeData);
        } else if ("-".equals(record.getAnnotations()[2])) {
          negRangesAccum.addRangeData(record.getSequenceName(), rangeData);
        } else {
          throw new IOException("BED record does not indicate strand as '-' or '+' in column 6: " + record);
        }
      }
    }
    mPosRanges = posRangesAccum.getReferenceRanges();
    mNegRanges = negRangesAccum.getReferenceRanges();
  }

  private static ProbeCounter checkList(RangeList<ProbeCounter> list, SAMRecord record, SAMRecord mate, int tolerance, PositionAndStrandChecker checker) {
    if (list == null) {
      return null;
    }
    if (!checker.checkStrand(record)) {
      return null;
    }
    int index = checker.getStartDataIndex(record, list);
    RangeList.RangeView<ProbeCounter> data = list.getFullRangeList().get(index);
    while (recordOverlap(record, data, tolerance)) {
      if (data != null && data.hasRanges()) {
        final ProbeCounter counter = data.getEnclosingRanges().get(0).getMeta();
        if (checker.checkPosition(record, data)) {
          checker.stripRecord(record, mate, data);
          counter.increment();
          return counter;
        }
      }
      if (list.getFullRangeList().size() > ++index) {
        data = list.getFullRangeList().get(index);
      } else {
        data = null;
      }
    }
    return null;
  }

  private static <T> boolean recordOverlap(SAMRecord rec, RangeList.RangeView<T> data, int tolerance) {
    if (data == null) {
      return false;
    }
    final int dataMin; // 0-based inclusive start
    if (data.getStart() - tolerance > data.getStart()) {
      dataMin = Integer.MIN_VALUE;
    } else {
      dataMin = data.getStart() - tolerance;
    }
    final int dataMax; // 0-based exclusive end
    if (data.getEnd() + tolerance < data.getEnd()) {
      dataMax = Integer.MAX_VALUE;
    } else {
      dataMax = data.getEnd() + tolerance;
    }
    final int alignmentStart = rec.getAlignmentStart() - 1;
    if (alignmentStart >= dataMin && alignmentStart < dataMax) {
      return true;
    }
    final int alignmentEnd = rec.getAlignmentEnd();
    if (alignmentEnd > dataMin && alignmentEnd <= dataMax) {
      return true;
    }
    return false;
  }

  private static final class ReadStatus {
    final String mName;
    boolean mStripped = false;
    boolean mNegative = false;
    private ReadStatus(String name, boolean stripped, boolean negative) {
      mName = name;
      setStatus(stripped, negative);
    }
    public int hashCode() {
      return mName.hashCode();
    }
    public boolean equals(Object that) {
      return that instanceof ReadStatus && mName.equals(((ReadStatus) that).mName);
    }
    private boolean stripped() {
      return mStripped;
    }
    private boolean negative() {
      return mNegative;
    }
    private void setStatus(boolean stripped, boolean negative) {
      mStripped = stripped;
      mNegative = negative;
    }
  }
}
