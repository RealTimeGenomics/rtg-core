/*
 * Copyright (c) 2014. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */
package com.rtg.simulation;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintStream;
import java.io.Writer;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.LoggedCli;
import com.rtg.reader.IndexFile;
import com.rtg.reader.PrereadNames;
import com.rtg.reader.PrereadNamesInterface;
import com.rtg.reader.ReaderUtils;
import com.rtg.reader.SdfId;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.reader.SourceTemplateReadWriter;
import com.rtg.sam.DefaultSamFilter;
import com.rtg.sam.SamFilterOptions;
import com.rtg.sam.SamFilterParams;
import com.rtg.sam.SamUtils;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.format.FormatReal;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.LineWriter;
import com.rtg.util.io.LogStream;

import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloseableIterator;

/**
 * Class to calculate
 *
 * <ul>
 * <li><b>correct</b> = read mapped to position from which it came,
 * accommodating known indels</li>
 * <li><b>incorrect</b> = read mapped uniquely to a different location with a
 * poor alignment score, not mapped, read mapped to multiple locations with same
 * alignment score <code>(topequal)</code> that does not include the original
 * location</li>
 * <li><b>better</b> = maps to one or more locations with a higher alignment
 * score (mapped with lower error than was introduced in the simulated data)</li>
 * </ul>
 *
 * Calculates accuracy / sensitivity by following formula <br>
 * <ul>
 * <li><b> Accuracy (Precision)</b> - correct / total mapped reads</li>
 * <li><b> Sensitivity (Recall)</b> - correct / total reads</li>
 * </ul>
 *
 * Calculates <b>lenient</b> accuracy / sensitivity by following formula <br>
 * <ul>
 * <li><b> Accuracy (Precision)</b> - (correct + better) / total mapped reads</li>
 * <li><b> Sensitivity (Recall)</b> - (correct + better) / total reads</li>
 * </ul>
 *
 */
public class ReadMappingAccuracy extends LoggedCli {

  private static final String MODULE_NAME = "readsimeval";

  private static final int MAX_WARNINGS = 10;

  private int mReadIdCounter = 0;
  private final HashMap<Integer, Integer> mNameMap = new HashMap<>();
  private SimulatedReadNameParser mParser = null;
  private SdfId[] mTemplateMap = null;
  private SdfId mOriginalReference = null;
  private PrereadNamesInterface mLeftNames = null;
  private PrereadNamesInterface mRightNames = null;
  private PrereadNamesInterface mLeftSuffixes = null;
  private PrereadNamesInterface mRightSuffixes = null;
  private ReadMappingAccuracyParams mParams = null;

  ReadMappingAccuracyReadStats mLeftStats;
  private ReadMappingAccuracyReadStats mRightStats;
  private final ReadMappingRoc mRocAs = new ReadMappingRoc();
  private final ReadMappingRoc mRocGs = new ReadMappingRoc();
  private final ReadMappingRoc mRocMapq = new ReadMappingRoc(false);

  private int mChimeras = 0;
  private int mDupes = 0;

  @Override
  protected void initFlags() {
    ReadMappingAccuracyParams.initFlags(mFlags);
  }

  @Override
  protected File outputDirectory() {
    return (File) mFlags.getValue(CommonFlags.OUTPUT_FLAG);
  }

  @Override
  protected int mainExec(OutputStream out, LogStream log) throws IOException {
    int exitCode = 0;
    final PrintStream psOut = new PrintStream(out);
    try {
      process(psOut);
    } catch (final InvalidParamsException e) {
      //treat this as an error in the arguments passed to the process
      mFlags.error(mFlags.getInvalidFlagMsg());
      cleanDirectory();
      return 1;
    } catch (final IllegalArgumentException ex) {
      Diagnostic.error(ex.getMessage());
      exitCode = 1;
    } finally {
      psOut.flush();
    }
    return exitCode;
  }

  private void setParser(String readName) throws IOException {
    if (mParser == null) {
      if (mParams.mutationsVcf() != null) {
        final MutatedSampleOffsets offsets = MutatedSampleOffsets.getOffsets(mParams.mutationsVcf(), mParams.sample());
        if (offsets == null) {
          throw new NoTalkbackSlimException("Overlapping variants in mutations VCF file.");
        }
        mParser = SimulatedReadNameParserFactory.getParsers(readName, offsets)[0];
      } else {
        mParser = SimulatedReadNameParserFactory.getParser(readName);
      }
    }
  }

  private int mIgnoredAbsentTemplate = 0;
  private boolean mHaveWarnedNoTemplateInSam = false;
  private boolean mHaveWarnedTemplateMismatchInSam = false;

  private void warnNoTemplateInSam() {
    if (!mHaveWarnedNoTemplateInSam) {
      Diagnostic.warning("Missing template SDF ID in SAM header. Cannot verify evaluation is against correct template.");
      mHaveWarnedNoTemplateInSam = true;
    }
  }

  private void warnTemplateMismatchInSam() {
    if (!mHaveWarnedTemplateMismatchInSam) {
      Diagnostic.warning("Encountered mappings against a template for which we cannot determine correct location.");
      mHaveWarnedTemplateMismatchInSam = true;
    }
    mIgnoredAbsentTemplate++;
  }

  private static final Pattern SOFT_CLIP_CIGAR_START_PATTERN = Pattern.compile("^([0-9]+)S([0-9]+).*");
  private SoftClipCigarParser mSoftClipParser = null;
  private int mSoftClipParseWarnings = 0;

  int determineAlignmentStart(SAMRecord rec, long originalLocation) {
    final Matcher softClipMatcher = SOFT_CLIP_CIGAR_START_PATTERN.matcher(rec.getCigarString());
    if (softClipMatcher.matches()) {
      //this read starts with a soft clip, so we now have to adjust the start position based on that.
      try {
        final int softClipLength = Integer.parseInt(softClipMatcher.group(1));

        //only bother inspecting cigars if the original location is "nearby" the aligned startposition (based on length of soft clipping)
        if (originalLocation > rec.getAlignmentStart() - softClipLength * 2 && originalLocation < rec.getAlignmentStart()) {
          int indelOffset = 0;
          if (mParser instanceof NewReadNameParser) {
            if (mSoftClipParser == null) {
              mSoftClipParser = new SoftClipCigarParser();
            }
            mSoftClipParser.setState(originalLocation, rec.getAlignmentStart(), ((NewReadNameParser) mParser).cigar(), rec.getCigarString());
            indelOffset = mSoftClipParser.parse();
          }
          return rec.getAlignmentStart() - softClipLength + indelOffset;
        } //else, not even in the ballpark. fall through.
      } catch (final NumberFormatException nfe) {
        if (mSoftClipParseWarnings++ < MAX_WARNINGS) {
          Diagnostic.warning("Unable to parse soft clip length in cigar: " + rec.getCigarString() + " of read " + rec.getReadName());
          if (mSoftClipParseWarnings == MAX_WARNINGS) {
            Diagnostic.warning("Further messages about soft clip length parsing errors will be suppressed.");
          }
        }
      }
    }
    return rec.getAlignmentStart();
  }

  // Cache the guid for each header to prevent having to reparse the headers for every single record.
  private final HashMap<SAMFileHeader, SdfId> mTemplateIds = new HashMap<>();
  SdfId getTemplateGuid(SAMFileHeader header) {
    if (!mTemplateIds.containsKey(header)) {
      mTemplateIds.put(header, SamUtils.getReferenceGuid(header));
    }
    return mTemplateIds.get(header);
  }

  void processRecord(SAMRecord rec, LineWriter recordsOut) throws IOException {
    final String readName = getReadName(rec);
    final int readLength = getReadLength(rec);
    setParser(readName);
    if (mParser == null || !mParser.setReadInfo(readName, readLength)) {
      throw new NoTalkbackSlimException("Read " + readName + " does not appear to be a simulated read");
    }
    if (mParser.isChimera()) {
//      System.err.println("Chimera mapped: " + rec.getSAMString());
      mChimeras++;
      return;
    } else if (mParser.isDuplicate()) {
//      System.err.println("Dupe mapped: " + rec.getSAMString());
      mDupes++;
    }
    final SdfId mappingRefGuid = getTemplateGuid(rec.getHeader());
    if (!mappingRefGuid.available()) {
      warnNoTemplateInSam();
    } else if ((mOriginalReference != null) && (!mappingRefGuid.check(mOriginalReference))
        || (mTemplateMap != null) && (!mappingRefGuid.check(mTemplateMap[mParser.templateSet()]))) {
      warnTemplateMismatchInSam();
      return;
    }
    final int id = getId(mParser.readId());
    final boolean isFirst = !mParams.isPaired() || rec.getFirstOfPairFlag();
    final ReadMappingAccuracyReadStats currentStats = isFirst ? mLeftStats : mRightStats;

    if (mParams.isPaired() && rec.getProperPairFlag()) {
      currentStats.mated(id);
    }
    final boolean isMultiple = currentStats.isMapped(id) || checkIfMultiple(rec); // Have we already seen it ourselves or did they fess-up to being multiply mapped
    currentStats.mapped(id);
    final double weight = getWeight(rec);
    if (isMultiple) {
      currentStats.multiple(id);
    }
    final String currentTemplateName = rec.getReferenceName();
    final long originalLocation = mParser.templatePosition();
    final long currentLocation = determineAlignmentStart(rec, originalLocation);
    final Integer as = rec.getIntegerAttribute(SamUtils.ATTRIBUTE_ALIGNMENT_SCORE);
    final int generatedScore = generatedScore(mParams.misMatchPenalty(), mParams.gapOpeningPenalty());
    if (currentTemplateName.equals(mParser.templateName()) && Math.abs(currentLocation - originalLocation) <= mParams.variance()) {
      // This record mapped the read to the correct location
      currentStats.found(id);
      if (as != null && mParams.scoreHistograms()) {
        mRocAs.addTp(as, weight);
        mRocGs.addTp(generatedScore, weight);
      }
      if (mParams.mapQHistogram() || mParams.mapQRoc()) {
        mRocMapq.addTp(rec.getMappingQuality(), weight);
      }
      if (mParams.dumpRecords()) {
        recordsOut.writeln(rec.getSAMString().trim() + "\t" + SamUtils.ATTRIBUTE_READ_ACCURACY_STATUS + ":Z:Correct");
      }
    } else { // The mapping was incorrect
      final Integer nm = rec.getIntegerAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES);
      if (as != null) {
        if (mParams.scoreHistograms()) {
          mRocAs.addFp(as, weight);
          mRocGs.addFp(generatedScore, weight);
        }
        if (as < generatedScore) {
          currentStats.better(id);
        }
      } else if (nm != null) {
        if (nm < mParser.numMismatches()) {
          currentStats.better(id);
        }
      }
      if (mParams.mapQHistogram() || mParams.mapQRoc()) {
        mRocMapq.addFp(rec.getMappingQuality(), weight);
      }
      if (mParams.dumpRecords()) {
        recordsOut.writeln(rec.getSAMString().trim() + "\t" + SamUtils.ATTRIBUTE_READ_ACCURACY_STATUS + ":Z:Incorrect");
      }
    }
  }

  void init() throws IOException, InvalidParamsException {
    mParams = new ReadMappingAccuracyParams(mFlags);
    mLeftStats = new ReadMappingAccuracyReadStats(mParams.isPaired() ? getNumReads(ReaderUtils.getLeftEnd(mParams.readDirectory())) : getNumReads(mParams.readDirectory()));
    mRightStats = mParams.isPaired() ? new ReadMappingAccuracyReadStats(getNumReads(ReaderUtils.getRightEnd(mParams.readDirectory()))) : null;
    if (mParams.mutationsVcf() != null) {
      mOriginalReference = SourceTemplateReadWriter.readMutationMap(mParams.isPaired() ? ReaderUtils.getLeftEnd(mParams.readDirectory()) : mParams.readDirectory());
      if (mOriginalReference == null) {
        Diagnostic.warning("No original reference map present in reads SDF. Cannot verify evaluation is against correct reference.");
      } else {
        Diagnostic.userLog("Evaluating accuracy of mapping reads drawn from mutated reference");
      }
    } else {
      mTemplateMap = SourceTemplateReadWriter.readTemplateMap(mParams.isPaired() ? ReaderUtils.getLeftEnd(mParams.readDirectory()) : mParams.readDirectory());
      if (mTemplateMap == null) {
        Diagnostic.warning("No template map present in reads SDF. Cannot verify evaluation is against correct template.");
      } else {
        Diagnostic.userLog("Evaluating accuracy of mapping reads drawn directly from reference");
      }
    }
  }

  private void process(final PrintStream out) throws IOException, InvalidParamsException {
    long totalRecords = 0;
    long skippedRecords = 0;
    long unmappedRecords = 0;
    boolean doneHeader = false;
    init();
    try (LineWriter recordsOut = mParams.dumpRecords() ? new LineWriter(new OutputStreamWriter(FileUtils.createOutputStream(new File(outputDirectory(), "mappings.sam")))) : null) {
      final SamFilterParams filterParams = SamFilterOptions.makeFilterParamsBuilder(mFlags).excludeUnmapped(true).excludeUnplaced(true).create();

      for (final File f : mParams.samFiles()) {
        try (final SamReader reader = new SAMFileReader(FileUtils.createInputStream(f, false))) {
          if (recordsOut != null && !doneHeader) {
            recordsOut.write(reader.getFileHeader().getTextHeader());
            doneHeader = true;
          }
          final CloseableIterator<SAMRecord> itr = reader.iterator();
          try {
            while (itr.hasNext()) {
              totalRecords++;
              if (totalRecords % 1000000 == 0) {
                Diagnostic.progress("Processed " + totalRecords + " SAM records");
              }
              final SAMRecord rec = itr.next();
              if (rec.getReadUnmappedFlag()) {
                unmappedRecords++;
              } else if (!DefaultSamFilter.acceptRecord(filterParams, rec)) {
                // ignore records as specified
                skippedRecords++;
              } else {
                processRecord(rec, recordsOut);
              }
            }
          } finally {
            itr.close();
          }
        }
      }
    }
    Diagnostic.progress("Processed " + totalRecords + " SAM records in total");
    printResults(out, totalRecords, unmappedRecords, skippedRecords);
  }

  int generatedScore(int mismatchPenalty, int gapOpeningPenalty) {
    return mParser.substitutions() * mismatchPenalty + (mParser.insertions() * gapOpeningPenalty) + (mParser.deletions() * gapOpeningPenalty);
  }

  // Get internal ID from external ID (which may be in a much larger space)
  private int getId(int readName) {
    // this was done to support partial read mappings
    if (mNameMap.containsKey(readName)) {
      return mNameMap.get(readName);
    } else {
      mNameMap.put(readName, mReadIdCounter);
      return mReadIdCounter++; // increment counter before returning
    }
  }

  private void printResults(PrintStream out, long totalRecords, long unmappedRecords, long skippedRecords) throws IOException {
    final File summary = new File(outputDirectory(), CommonFlags.SUMMARY_FILE);
    try (LineWriter summaryout = new LineWriter(new OutputStreamWriter(FileUtils.createOutputStream(summary)))) {
      final FormatReal formatter = new FormatReal(3, 2);
      summaryout.writeln("Total SAM records = " + totalRecords);
      summaryout.writeln("Unmapped records = " + unmappedRecords);
      summaryout.writeln("Filtered records = " + skippedRecords);
      if (mParams.isPaired()) {
        printPairedEndStats(summaryout, formatter);
      } else {
        printSingleEndStats(summaryout, formatter);
      }
      if (mChimeras > 0) {
        summaryout.writeln("Chimeras (both sides) = " + mChimeras);
      }
      if (mDupes > 0) {
        summaryout.writeln("Dupes (both sides) = " + mDupes);
      }
      printCommonStats();
    }
    out.print(FileUtils.fileToString(summary));
  }

  private void printCommonStats() throws IOException {
    if (mParams.scoreHistograms()) {
      try (Writer scoreout = new BufferedWriter(new OutputStreamWriter(FileUtils.createOutputStream(new File(outputDirectory(), "score_hist.tsv"))))) {
        scoreout.write("Alignment Score Distribution" + StringUtils.LS);
        scoreout.write(mRocAs.getDistribution());
        if (mParams.verbose()) {
          scoreout.write(StringUtils.LS);
          scoreout.write("Generated Score Distribution" + StringUtils.LS);
          scoreout.write(mRocGs.getDistribution());
        }
      }
    }
    if (mParams.mapQHistogram()) {
      try (Writer mapqout = new BufferedWriter(new OutputStreamWriter(FileUtils.createOutputStream(new File(outputDirectory(), "mapq_hist.tsv"))))) {
        mapqout.write("MAPQ Score Distribution" + StringUtils.LS);
        mapqout.write(mRocMapq.getDistribution());
      }
    }
    if (mParams.mapQRoc()) {
      try (Writer rocout = new BufferedWriter(new OutputStreamWriter(FileUtils.createOutputStream(new File(outputDirectory(), "mapq_roc.tsv"))))) {
        rocout.write("#MAPQ ROC" + StringUtils.LS);
        rocout.write(mRocMapq.getRoc(mLeftStats.length() * (mParams.isPaired() ? 2 : 1)));
      }
    }
  }

  private static class SummaryStats2 {
    private final String mLabel;
    long mTotal = 0;
    long mAndUnique = 0;
    long mAndMultiple = 0;

    SummaryStats2(String label) {
      mLabel = label;
    }
  }

  private static class SummaryStats {
    private final String mLabel;
    private final SummaryStats2 mCorrect = new SummaryStats2("correctly");
    private final SummaryStats2 mIncorrect = new SummaryStats2("incorrectly");
    private final SummaryStats2 mBetter = new SummaryStats2("incorrectly with lower alignment score");
    SummaryStats(String label) {
      mLabel = label;
    }
  }

  private static void updateStats(ReadMappingAccuracyReadStats readStats, int i, SummaryStats2 stats) {
    stats.mTotal++;
    if (readStats.isMultiple(i)) {
      stats.mAndMultiple++;
    } else {
      stats.mAndUnique++;
    }
  }

  private static void updateStats(ReadMappingAccuracyReadStats readStats, int i, SummaryStats stats) {
    if (readStats.isFound(i)) {
      updateStats(readStats, i, stats.mCorrect);
    } else if (readStats.isMapped(i)) {
      updateStats(readStats, i, stats.mIncorrect);
      if (readStats.isBetter(i)) {
        updateStats(readStats, i, stats.mBetter);
      }
    } // Else unmapped
  }

  private void printSingleEndStats(LineWriter summaryout, final FormatReal formatter) throws IOException {
    final long totalReads = mLeftStats.length();
    final SummaryStats stats = new SummaryStats("reads");
    long mappedReads = 0;
    long unmappedReads = 0;

    for (int i = 0; i < mLeftStats.length(); i++) {
      updateStats(mLeftStats, i, stats);
      if (mLeftStats.isMapped(i)) {
        mappedReads++;
      } else {
        unmappedReads++;
      }
    }
    if (mIgnoredAbsentTemplate > 0) {
      summaryout.writeln("Records ignored due to simulation/mapping template mismatch = " + mIgnoredAbsentTemplate);
    }
    summaryout.writeln("Total reads = " + totalReads);
    summaryout.writeln("mapped reads = " + mappedReads);
    assert unmappedReads == (totalReads - mappedReads);
    summaryout.writeln("unmapped reads = " + unmappedReads);

    for (final SummaryStats2 stats2 : new SummaryStats2[] {stats.mCorrect, stats.mIncorrect, stats.mBetter}) {
      summaryout.writeln(stats.mLabel + " mapped " + stats2.mLabel + " = " + stats2.mTotal);
    }

    final double totalCorrect = stats.mCorrect.mTotal;
    final double totalIncorrect = stats.mIncorrect.mTotal; // Includes better
    final double totalBetter = stats.mBetter.mTotal;
    final double totalMapped = totalCorrect + totalIncorrect;
    printAccuracy(summaryout, formatter, totalReads, totalCorrect, totalBetter, totalMapped);

    if (mParams.verbose()) {
      for (final SummaryStats2 stats2 : new SummaryStats2[] {stats.mCorrect, stats.mIncorrect, stats.mBetter}) {
        summaryout.writeln(stats.mLabel + " mapped " + stats2.mLabel + " and unique = " + stats2.mAndUnique);
        summaryout.writeln(stats.mLabel + " mapped " + stats2.mLabel + " and multiple = " + stats2.mAndMultiple);
      }
    }
  }

  private void printPairedEndStats(LineWriter summaryout, final FormatReal formatter) throws IOException {
    final long total = mLeftStats.length();
    final SummaryStats leftStats = new SummaryStats("left reads");
    final SummaryStats rightStats = new SummaryStats("right reads");
    long mated = 0;
    long unmated = 0;

    for (int i = 0; i < mLeftStats.length(); i++) {
      updateStats(mLeftStats, i, leftStats);
      updateStats(mRightStats, i, rightStats);
      if (mLeftStats.isMapped(i) || mRightStats.isMapped(i)) {
        if (mRightStats.isMated(i)) {
          mated++;
        } else {
          unmated++;
        }
      }
    }
    if (mIgnoredAbsentTemplate > 0) {
      summaryout.writeln("Records ignored due to simulation/mapping template mismatch = " + mIgnoredAbsentTemplate);
    }
    summaryout.writeln("Total pairs = " + total);
    summaryout.writeln("mated pairs = " + mated);
    summaryout.writeln("unmated pairs = " + unmated);

    for (final SummaryStats stats : new SummaryStats[] {leftStats, rightStats}) {
      for (final SummaryStats2 stats2 : new SummaryStats2[] {stats.mCorrect, stats.mIncorrect, stats.mBetter}) {
        if (mParams.verbose() || stats2 != stats.mBetter) { // Only output "better" stuff if verbose, it's probably broken across different mapping tools
          summaryout.writeln(stats.mLabel + " mapped " + stats2.mLabel + " = " + stats2.mTotal);
        }
      }
    }

    final double totalCorrect = leftStats.mCorrect.mTotal + rightStats.mCorrect.mTotal;
    final double totalIncorrect = leftStats.mIncorrect.mTotal + rightStats.mIncorrect.mTotal;
    final double totalBetter = leftStats.mBetter.mTotal + rightStats.mBetter.mTotal;
    final double totalMapped = totalCorrect + totalIncorrect; // better is included in incorrect
    printAccuracy(summaryout, formatter, total * 2, totalCorrect, totalBetter, totalMapped);

    if (mParams.verbose()) {
      for (final SummaryStats stats : new SummaryStats[] {leftStats, rightStats}) {
        for (final SummaryStats2 stats2 : new SummaryStats2[] {stats.mCorrect, stats.mIncorrect, stats.mBetter}) {
          if (mParams.verbose() || stats2 != stats.mBetter) { // Only output "better" stuff if verbose, it's probably broken across different mapping tools
            summaryout.writeln(stats.mLabel + " mapped " + stats2.mLabel + " and unique = " + stats2.mAndUnique);
            summaryout.writeln(stats.mLabel + " mapped " + stats2.mLabel + " and multiple = " + stats2.mAndMultiple);
          }
        }
      }
    }
  }

  private void printAccuracy(LineWriter summaryout, final FormatReal formatter, final long total, final double totalCorrect, final double totalBetter, final double totalMapped) throws IOException {
    if (mParams.verbose()) {
      final double lenientAccuracy = (totalCorrect + totalBetter) / totalMapped * 100;
      final double lenientSensitivity = (totalCorrect + totalBetter) / total * 100;
      summaryout.writeln("Lenient Accuracy (Precision) = (correct + mapped incorrectly with lower alignment score) / total mapped reads = " + formatter.format(lenientAccuracy));
      summaryout.writeln("Lenient Sensitivity (Recall) = (correct + mapped incorrectly with lower alignment score) / total reads = " + formatter.format(lenientSensitivity));
    }
    final double accuracy = 100.0 * totalCorrect / totalMapped;
    final double sensitivity = 100.0 * totalCorrect / total;
    summaryout.writeln("Accuracy (Precision) = correct / total mapped reads = " + formatter.format(accuracy));
    summaryout.writeln("Sensitivity (Recall) = correct / total reads = " + formatter.format(sensitivity));
  }

  // Should be 1 / number of mappings produced for the read
  private double getWeight(SAMRecord rec) {
    // this is for RTG, maybe TMAP
    final Integer res = rec.getIntegerAttribute(SamUtils.ATTRIBUTE_IH);
    if (res != null) {
      return 1.0 / res;
    }

    // While X0 indicates the number of best scoring alignments - they
    // are not necessarily all output to the SAM file. Allowing a reduction
    // in weight for these cases will unfairly benefit BWA during the analysis,
    // as mappings more likely to be wrong will be given less total weight. We can
    // only fairly distribute the weight among records actually output.
    // this is for BWA
    final Integer x0 = rec.getIntegerAttribute(SamUtils.ATTRIBUTE_BWA_NUM_BEST_HITS);
    if (x0 != null) {
      return 1.0 / x0;
    }

    // Otherwise assume a single mapping
    return 1.0;
  }

  private boolean checkIfMultiple(SAMRecord rec) {
    // this is for RTG, maybe TMAP
    final Integer res = SamUtils.getNHOrIH(rec);
    if (res != null) {
      return res > 1;
    }

    // this is for BWA
    try {
      final Object o = rec.getCharacterAttribute(SamUtils.ATTRIBUTE_BWA_TYPE);
      if (o != null) {
        final Character xt = (Character) o;
        if (Character.toLowerCase(xt) == 'u') {
          return false;
        }
      }
    } catch (SAMException se) {
      //if it's the wrong type, this may not be BWA - just ignore for now.
    }

    final Integer x0 = rec.getIntegerAttribute(SamUtils.ATTRIBUTE_BWA_NUM_BEST_HITS);
    return x0 != null && x0 > 1;
  }

  private int getReadLength(SAMRecord rec) {
    return rec.getCigar().getReferenceLength();
  }

  String leftName(long readNum) throws IOException {
    if (mLeftNames == null) {
      final File preread = mParams.isPaired() ? ReaderUtils.getLeftEnd(mParams.readDirectory()) : mParams.readDirectory();
      mLeftNames = new PrereadNames(preread, LongRange.NONE);
      final IndexFile index = new IndexFile(preread);
      if (index.hasSequenceNameSuffixes()) {
        mLeftSuffixes = new PrereadNames(preread, LongRange.NONE, true);
      }
    }
    final String name;
    if (mLeftSuffixes != null) {
      name =  mLeftNames.name(readNum) + " " + mLeftSuffixes.name(readNum).trim();
    } else {
      name =  mLeftNames.name(readNum);
    }
    return name;
  }

  String rightName(long readNum) throws IOException {
    if (mRightNames == null) {
      final File preread = mParams.isPaired() ? ReaderUtils.getRightEnd(mParams.readDirectory()) : mParams.readDirectory();
      mRightNames = new PrereadNames(preread, LongRange.NONE);
      final IndexFile index = new IndexFile(preread);
      if (index.hasSequenceNameSuffixes()) {
        mRightSuffixes = new PrereadNames(preread, LongRange.NONE, true);
      }
    }
    final String name;
    if (mRightSuffixes != null) {
      name =  mRightNames.name(readNum) + " " + mRightSuffixes.name(readNum).trim();
    } else {
      name =  mRightNames.name(readNum);
    }
    return name;
  }

  private String getReadName(SAMRecord rec) throws IOException {
    try {
      final long readNum = Long.parseLong(rec.getReadName());
      if (mParams.isPaired()) {
        if (rec.getFirstOfPairFlag()) {
          return leftName(readNum);
        } else {
          return rightName(readNum);
        }
      } else {
        return leftName(readNum);
      }
    } catch (final NumberFormatException nfe) {
      return rec.getReadName();
    }
  }

  @Override
  public String moduleName() {
    return MODULE_NAME;
  }

  /** @param args command line arguments */
  public static void main(String[] args) {
    new ReadMappingAccuracy().mainExit(args);
  }

  private int getNumReads(File seqDir) throws IOException {
    final int result;
    try (SequencesReader dsr = SequencesReaderFactory.createDefaultSequencesReader(seqDir)) {
      result = (int) dsr.numberSequences();
    }
    return result;
  }
}
