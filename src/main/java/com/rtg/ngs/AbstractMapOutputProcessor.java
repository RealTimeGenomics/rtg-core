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
package com.rtg.ngs;

import static com.rtg.index.hash.ngs.ReadDecoder.PAIRED_END;
import static com.rtg.ngs.ReadStatusTracker.SHORT_FIRST;
import static com.rtg.ngs.ReadStatusTracker.SHORT_SECOND;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.calibrate.SamCalibrationInputs;
import com.rtg.index.hash.ngs.OutputProcessor;
import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.HashingRegion;
import com.rtg.launcher.ISequenceParams;
import com.rtg.launcher.globals.CoreGlobalFlags;
import com.rtg.launcher.globals.GlobalFlags;
import com.rtg.ngs.blocking.MapQScoringReadBlocker;
import com.rtg.ngs.blocking.MapQScoringReadBlockerSynch;
import com.rtg.ngs.tempstage.AbstractTempFileWriter;
import com.rtg.ngs.tempstage.PairedTempFileWriterImpl;
import com.rtg.ngs.tempstage.SingleEndTempFileWriter;
import com.rtg.reader.NamesInterface;
import com.rtg.reader.SequencesReader;
import com.rtg.sam.SamFilterParams;
import com.rtg.sam.SamMerger;
import com.rtg.sam.SamUtils;
import com.rtg.util.IORunnable;
import com.rtg.util.NullStreamUtils;
import com.rtg.util.SimpleThreadPool;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.intervals.Range;
import com.rtg.util.intervals.RangeList;
import com.rtg.util.intervals.RangeMeta;
import com.rtg.util.intervals.SimpleRangeMeta;
import com.rtg.util.io.FileUtils;
import com.rtg.variant.sv.ReadGroupStatsCalculator;
import com.rtg.variant.sv.UnmatedAugmenter;

import htsjdk.samtools.SAMRecord;

/**
 * Class to provide common features between single and paired end output processors for map command.
 * See <a href="http://intranet.nz.realtimegenomics.com/wiki/index.php?n=RTG.Mapping"><code>wiki</code>:RTG.Mapping</a> for more information about workflow of this class and its subclasses.
 *
 */
@TestClass(value = {"com.rtg.ngs.SamSingleEndOutputProcessorTest", "com.rtg.ngs.TopNPairedEndOutputProcessorSyncTest", "com.rtg.ngs.DummyMapOutputProcessorTest"})
public abstract class AbstractMapOutputProcessor implements OutputProcessor {
  protected static final String TEMP_SAM_ALIGNMENT_NAME = "TEMP_SAM_";
  private static final String TEMP_SAM_UNMATED_ALIGNMENT_NAME = "TEMP_SAM_UNMATED_";
  protected final SharedResources mSharedResources;
  protected final NgsParams mParams;
  protected final UnmatedAugmenter.Merger mAugmenterMerger;
  protected final ReadGroupStatsCalculator.Merger mStatsMerger;
  protected final MapReportData.Merger mReportMerger;
  protected final ReadStatusTracker mUnmappedTracker;
  protected final List<HashingRegion> mRegions;
  private final boolean mPaired;

  /**
   * @param param the parameters
   * @param stats statistics container
   * @param paired whether in paired or single end mode
   * @param outputUnmated if we should output unmated. Ignored if {@code paired} is false
   * @throws IOException if an IO error occurs
   */
  public AbstractMapOutputProcessor(NgsParams param, MapStatistics stats, boolean paired, boolean outputUnmated) throws IOException {
    final int sequences = (int) param.buildFirstParams().numberSequences();
    mSharedResources = SharedResources.generateSharedResources(param);
    mReportMerger = new MapReportData.Merger();
    mParams = param;
    mUnmappedTracker = new ReadStatusTrackerSync(sequences, stats);
    if (param.outputParams().ignoreShort()) {
      setShortReadFlags(mUnmappedTracker, param, paired);
    }
    if (paired && outputUnmated) {
      mAugmenterMerger = mParams.outputParams().svprep() ? new UnmatedAugmenter.Merger() : null;
      mStatsMerger = mParams.outputParams().svprep() ? new ReadGroupStatsCalculator.Merger() : null;
    } else {
      mAugmenterMerger = null;
      mStatsMerger = null;
    }
    mPaired = paired;
    mRegions = new ArrayList<>();
  }

  // Record which reads are shorter than the indexable length - these may be ignored from statistics and output
  private static void setShortReadFlags(ReadStatusTracker tracker, NgsParams param, boolean paired) throws IOException {
    final int threshold = minReadLength(param.maskParams());
    Diagnostic.developerLog("Flagging all reads shorter than " + threshold);
    int l = 0;
    int r = 0;
    final SequencesReader r1 = param.buildFirstParams().reader();
    final int len1 = (int) r1.numberSequences();
    for (int i = 0; i < len1; i++) {
      if (r1.length(i) < threshold) {
        tracker.addStatus(i, SHORT_FIRST);
        l++;
      }
    }
    if (paired) {
      final SequencesReader r2 = param.buildSecondParams().reader();
      final int len2 = (int) r2.numberSequences();
      for (int i = 0; i < len2; i++) {
        if (r2.length(i) < threshold) {
          tracker.addStatus(i, SHORT_SECOND);
          r++;
        }
      }
    }
    Diagnostic.developerLog("Flagged all reads shorter than " + threshold + "nt. " + l + " left, " + r + " right");
  }
  private static int minReadLength(NgsMaskParams mp) {
    try {
      return mp.getWordSize();
    } catch (UnsupportedOperationException e) {
      return 1; // For any unknown word size options, only ignore zero length reads
    }
  }

  protected abstract FilterConcatIntermediateFiles filterConcatNonMated(MapQScoringReadBlocker blockerLeft, MapQScoringReadBlocker blockerRight, File[] tempFiles, SingleEndTopRandomImplementation.HitRecord[] hitsToKeep, NamesInterface templateNames, File outFile) throws IOException;

  private static final class FileAndStream {
    private final File mFile;
    private final OutputStream mStream;

    private FileAndStream(File file, OutputStream stream) {
      mFile = file;
      mStream = stream;
    }

  }

  protected FilterConcatIntermediateFiles writeUnmapped(boolean isFinalUnmapped, boolean suppressSam, boolean unfiltered) throws IOException {
    mUnmappedTracker.preProcessUnMappedStatistics(mPaired);
    final ArrayList<File> retOutput = new ArrayList<>();
    final ArrayList<File> retCalibrate = new ArrayList<>();
    final FileAndStream unmappedOut = getUnmappedFileAndStream(isFinalUnmapped, suppressSam, "unrescued");
    Diagnostic.progress("UnmappedProcess: Starting 1 Jobs");
    try (final OutputStream unmappedOutStream = unmappedOut.mStream) {
      final UnmappedSamRecordFactory unmappedSamRecordFactory = new UnmappedSamRecordFactory(mParams, mSharedResources);
      try (UnmappedSamAlignmentWriter unmappedNoPositionSamWriter = new UnmappedSamAlignmentWriter(mParams.outputParams().tempFilesDirectory(), mSharedResources.getHeader())) {
        if (!unfiltered) {
          if (mPaired) {
            unmappedSamRecordFactory.setAugmenterInfo(mAugmenterMerger, mStatsMerger);
          }
          unmappedSamRecordFactory.setMapReportData(mReportMerger);
        }

        final boolean bam = mParams.outputParams().bam();
        unmappedNoPositionSamWriter.initialiseUnmapped(unmappedOutStream, bam, mAugmenterMerger == null, isFinalUnmapped);
        if (!isFinalUnmapped) {
          //alignment file unification mode
          final int numberIntermediateFiles = AbstractMulticoreFilterConcat.numberIntermediateFiles(mParams.numberThreads(), mRegions.size());
          final Range[] regionRanges = AbstractMulticoreFilterConcat.groupRegions(mRegions.size(), numberIntermediateFiles);
          final FileAndStream[] outputFiles = new FileAndStream[regionRanges.length];
          final UnmappedSamAlignmentWriter[] outputWriters = new UnmappedSamAlignmentWriter[regionRanges.length];
          try {
            for (int i = 0; i < outputFiles.length; ++i) {
              outputFiles[i] = getUnmappedFileAndStream(isFinalUnmapped, suppressSam, Integer.toString(i));
              outputWriters[i] = new UnmappedSamAlignmentWriter(mParams.outputParams().tempFilesDirectory(), mSharedResources.getHeader());
              outputWriters[i].initialiseUnmapped(outputFiles[i].mStream, bam, mAugmenterMerger == null, true);
              retOutput.add(outputFiles[i].mFile);
              retCalibrate.add(SamSingleEndOutputProcessor.calibrateUnmappedFile(mParams, outputFiles[i].mFile));
            }
            final Map<Long, RangeList<UnmappedSamAlignmentWriter>> lookup = getUnmappedWriterReferenceLookup(mParams.searchParams(), mRegions, regionRanges, outputWriters, unmappedNoPositionSamWriter);
            final RecordToWriter rtw = new LookupRecordToWriter(lookup);
            writeUnmappedRecords(unmappedSamRecordFactory, rtw);
          } finally {
            closeUnmappedWriters(outputWriters);
          }
        } else {
          //normal single file operation
          final RecordToWriter rtw = new SimpleRecordToWriter(unmappedNoPositionSamWriter);
          writeUnmappedRecords(unmappedSamRecordFactory, rtw);
        }
      }
    }
    final File indexFile;
    final File calibrateFile;
    if (unmappedOut.mFile != null) {
      indexFile = SamSingleEndOutputProcessor.indexSamFile(mParams, unmappedOut.mFile, isFinalUnmapped, mSharedResources.getHeader().getSequenceDictionary().size());
      calibrateFile = SamSingleEndOutputProcessor.calibrateUnmappedFile(mParams, unmappedOut.mFile);
    } else {
      indexFile = null;
      calibrateFile = null;
    }
    Diagnostic.progress("UnmappedProcess: 1/1 Jobs Finished");
    retOutput.add(unmappedOut.mFile);
    retCalibrate.add(calibrateFile);
    return new FilterConcatIntermediateFiles(retOutput.toArray(new File[0]), retCalibrate.toArray(new File[0]), !isFinalUnmapped || indexFile == null ? null : new File[] {indexFile});
  }

  private void writeUnmappedRecords(UnmappedSamRecordFactory unmappedSamRecordFactory, RecordToWriter rtw) {
    for (int read = 0; read < mUnmappedTracker.getNumReads(); ++read) {
      final ReadStatusTracker.UnmappedStatus status = mUnmappedTracker.getUnmappedStatus(read, mPaired);
      switch (status) {
        case SINGLE_END_UNMAPPED:
        case LEFT_UNMAPPED:
          writeUnmappedRecord(unmappedSamRecordFactory, rtw, read, true, false);
          break;
        case RIGHT_UNMAPPED:
          writeUnmappedRecord(unmappedSamRecordFactory, rtw, read, false, false);
          break;
        case BOTH_UNMAPPED:
          writeUnmappedRecord(unmappedSamRecordFactory, rtw, read, true, true);
          writeUnmappedRecord(unmappedSamRecordFactory, rtw, read, false, true);
          break;
        case MAPPED:
          break;
        default:
          throw new UnsupportedOperationException();
      }
    }
  }

  private static void closeUnmappedWriters(UnmappedSamAlignmentWriter[] outputWriters) {
    for (UnmappedSamAlignmentWriter writer : outputWriters) {
      if (writer != null) {
        writer.close();
      }
    }
  }

  protected static Map<Long, RangeList<UnmappedSamAlignmentWriter>> getUnmappedWriterReferenceLookup(ISequenceParams searchParams, List<HashingRegion> regions, Range[] regionRanges, UnmappedSamAlignmentWriter[] outputWriters, UnmappedSamAlignmentWriter unmappedNoPositionSamWriter) throws IOException {
    final Map<Long, List<RangeMeta<UnmappedSamAlignmentWriter>>> temp = new HashMap<>();
    final Map<Long, RangeList<UnmappedSamAlignmentWriter>> lookup = new HashMap<>();
    for (int rangeIndex = 0; rangeIndex < regionRanges.length; ++rangeIndex) {
      final Range range = regionRanges[rangeIndex];
      for (int regionIndex = range.getStart(); regionIndex < range.getEnd(); ++regionIndex) {
        final HashingRegion r = regions.get(regionIndex);
        final long regionStartId = r.getStart() != HashingRegion.MISSING ? r.getStart() : 0;
        final long regionEndId = r.getEnd() == HashingRegion.MISSING ? searchParams.numberSequences() : r.getExclusiveEndId();
        for (long seq = regionStartId; seq < regionEndId; ++seq) {
          final int startPos;
          if (seq != regionStartId) {
            startPos = 0;
          } else {
            startPos = r.getStartClipPosition() != HashingRegion.MISSING ? (int) r.getStartClipPosition() : 0;
          }
          final int endPos;
          if (seq != regionEndId - 1) {
            endPos = searchParams.reader().length(seq);
          } else {
            endPos = r.getEndClipPosition() != HashingRegion.MISSING ? (int) r.getEndClipPosition() : searchParams.reader().length(seq);
          }
          final RangeMeta<UnmappedSamAlignmentWriter> rr = new SimpleRangeMeta<>(startPos, endPos, outputWriters[rangeIndex]);
          temp.computeIfAbsent(seq, k -> new ArrayList<>()).add(rr);
        }
      }
    }
    for (Map.Entry<Long, List<RangeMeta<UnmappedSamAlignmentWriter>>> entry : temp.entrySet()) {
      lookup.put(entry.getKey(), new RangeList<>(entry.getValue()));
    }
    lookup.put(-1L, new RangeList<>(new SimpleRangeMeta<>(-1, Integer.MAX_VALUE, unmappedNoPositionSamWriter)));
    return lookup;
  }

  private void writeUnmappedRecord(UnmappedSamRecordFactory unmappedSamRecordFactory, RecordToWriter rtw, int read, boolean first, boolean mateUnmapped) {
    final boolean ignore = mParams.outputParams().ignoreShort() && mUnmappedTracker.getStatus(read, first ? SHORT_FIRST : SHORT_SECOND);
    if (!ignore) {
      final SAMRecord rec = unmappedSamRecordFactory.unmappedResult(read, first, mUnmappedTracker.getXCAttribute(read, first), mateUnmapped);
      final UnmappedSamAlignmentWriter unmappedSamAlignmentWriter = rtw.getWriter(rec);
      unmappedSamAlignmentWriter.unmappedRecord(rec);
    }
  }

  private interface RecordToWriter {
    UnmappedSamAlignmentWriter getWriter(SAMRecord r);
  }

  private static final class LookupRecordToWriter implements RecordToWriter {

    private final Map<Long, RangeList<UnmappedSamAlignmentWriter>> mLookup;
    private LookupRecordToWriter(Map<Long, RangeList<UnmappedSamAlignmentWriter>> lookup) {
      mLookup = lookup;
    }

    @Override
    public UnmappedSamAlignmentWriter getWriter(SAMRecord r) {
      final long referenceIndex = (long) r.getReferenceIndex();
      final int loc = r.getAlignmentStart() - 1;
      final RangeList<UnmappedSamAlignmentWriter> unmappedSamAlignmentWriterRangeList = mLookup.get(referenceIndex);
      final List<UnmappedSamAlignmentWriter> unmappedSamAlignmentWriters = unmappedSamAlignmentWriterRangeList.find(loc);
      if (unmappedSamAlignmentWriters == null) {
        throw new RuntimeException("No unmapped writer found for region: " + r.getReferenceName() + ":" + r.getAlignmentStart() + " for read: " + r.getReadName() + " " + r.getReadString());
      }
      return unmappedSamAlignmentWriters.get(0);
    }
  }
  private static final class SimpleRecordToWriter implements  RecordToWriter {
    private final UnmappedSamAlignmentWriter mWriter;

    private SimpleRecordToWriter(UnmappedSamAlignmentWriter writer) {
      mWriter = writer;
    }

    @Override
    public UnmappedSamAlignmentWriter getWriter(SAMRecord r) {
      return mWriter;
    }
  }


  private FileAndStream getUnmappedFileAndStream(boolean isFinalUnmapped, boolean suppressSam, String tag) throws IOException {
    final FileAndStream unmappedOut;
    if (mParams.outputParams().bam()) {
      final File unmappedOutFile;
      if (isFinalUnmapped) {
        unmappedOutFile = mParams.outputParams().file(NgsOutputParams.UNMAPPED_BAM_FILE_NAME);
      } else {
        unmappedOutFile = File.createTempFile("TEMP_UNMAPPED_" + tag + "_", ".bam", mParams.outputParams().tempFilesDirectory());
      }
      //we override picards bam constructor to allow us to provide the block gzip for the unmapped file
      unmappedOut = new FileAndStream(unmappedOutFile, FileUtils.createOutputStream(unmappedOutFile, true, true));
    } else if (!suppressSam) {
      final File unmappedOutFile;
      if (isFinalUnmapped) {
        unmappedOutFile = mParams.outputParams().isCompressOutput()
                ? mParams.outputParams().resultStreamHandler().file(NgsOutputParams.UNMAPPED_SAM_FILE_NAME + FileUtils.GZ_SUFFIX)
                : mParams.outputParams().resultStreamHandler().file(NgsOutputParams.UNMAPPED_SAM_FILE_NAME);
      } else {
        unmappedOutFile = File.createTempFile("TEMP_UNMAPPED_" + tag + "_", SamUtils.SAM_SUFFIX + (mParams.outputParams().isCompressOutput() ? FileUtils.GZ_SUFFIX : ""), mParams.outputParams().tempFilesDirectory());
      }
      unmappedOut = new FileAndStream(unmappedOutFile, FileUtils.createOutputStream(unmappedOutFile, mParams.outputParams().isCompressOutput(), mParams.outputParams().isCompressOutput()));
    } else {
      unmappedOut = new FileAndStream(null, NullStreamUtils.getNullOutputStream()); // TODO avoid creating SAMRecords entirely in this case
    }
    return unmappedOut;
  }

  protected void whizBangUnify(final FilterConcatIntermediateFiles unmappedAlignments, final FilterConcatIntermediateFiles... alignmentFiles) throws IOException {
    final boolean keepTempFiles = GlobalFlags.isSet(CoreGlobalFlags.MAP_KEEP_TEMPORARY_FILES);
    //if present the unmapped array contains unmapped, and will be 1 longer than the rest (containing the unmapped with no position entries)
    final SamMerger merge = new SamMerger(mParams.outputParams().outputIndex(), mParams.outputParams().isCompressOutput(), mParams.legacyCigars(), mParams.numberThreads(), SamFilterParams.builder().create(), false, !keepTempFiles);
    final int numberIntermediateFiles = AbstractMulticoreFilterConcat.numberIntermediateFiles(mRegions.size(), mParams.numberThreads());
    final SimpleThreadPool stp = new SimpleThreadPool(Math.min(numberIntermediateFiles, AbstractMulticoreFilterConcat.MAX_FILTERCONCAT_THREADS), "Region-Merge", true);
    final List<File> regionFiles = new ArrayList<>();
    try {
      for (int i = 0; i < numberIntermediateFiles; ++i) {
        final int index = i;
        final String fileNameSuffix = mParams.outputParams().bam() ? SamUtils.BAM_SUFFIX
                : (mParams.outputParams().isCompressOutput() ? SamUtils.SAM_SUFFIX + FileUtils.GZ_SUFFIX
                : SamUtils.SAM_SUFFIX);
        final File outputFile = File.createTempFile("TEMP_REGION_" + index + "_", fileNameSuffix, mParams.outputParams().tempFilesDirectory());
        regionFiles.add(outputFile);
        stp.execute(new IORunnable() {
          @Override
          public void run() throws IOException {
            final ArrayList<File> toMerge = new ArrayList<>();
            for (FilterConcatIntermediateFiles arr : alignmentFiles) {
              if (arr != null) {
                toMerge.add(arr.getAlignmentFiles().get(index));
              }
            }
            if (unmappedAlignments != null) {
              toMerge.add(unmappedAlignments.getAlignmentFiles().get(index));
            }
            final SamCalibrationInputs inputs = new SamCalibrationInputs(toMerge, true);
            merge.mergeSamFiles(inputs.getSamFiles(), inputs.getCalibrationFiles(), outputFile, null, null, mSharedResources.getHeader().clone(), index == 0, false);
          }
        });
      }
    } finally {
      stp.terminate();
    }
    if (unmappedAlignments != null) {
      regionFiles.add(unmappedAlignments.getAlignmentFiles().get(unmappedAlignments.getAlignmentFiles().size() - 1));
    }
    final File finalOutputFile = mParams.outputParams().file(
            mParams.outputParams().bam() ? NgsOutputParams.ALIGNMENTS_BAM_FILE_NAME
                    : (mParams.outputParams().isCompressOutput()
                    ? NgsOutputParams.ALIGNMENTS_SAM_FILE_NAME + FileUtils.GZ_SUFFIX
                    : NgsOutputParams.ALIGNMENTS_SAM_FILE_NAME)
    );
    final File[] indexFiles = new File[regionFiles.size()];
    final File[] calibrationFiles = new File[regionFiles.size()];
    final ArrayList<Long> dataFileSizes = new ArrayList<>();
    for (int i = 0; i < regionFiles.size(); ++i) {
      final File f = regionFiles.get(i);
      dataFileSizes.add(f.length());
      indexFiles[i] = AbstractMulticoreFilterConcat.indexFileName(f, mParams.outputParams().bam());
      calibrationFiles[i] = new File(f.getParent(), f.getName() + CommonFlags.RECALIBRATE_EXTENSION);
    }
    FileUtils.copyRaw(finalOutputFile, !keepTempFiles, regionFiles.toArray(new File[0]));
    if (mParams.outputParams().outputIndex() && (mParams.outputParams().isCompressOutput() || mParams.outputParams().bam())) {
      AbstractMulticoreFilterConcat.mergeIndexes(mParams, finalOutputFile, indexFiles, dataFileSizes);
    }
    if (mParams.outputParams().calibrate()) {
      AbstractMulticoreFilterConcat.mergeCalibration(finalOutputFile, mParams.outputParams().calibrateRegions(), calibrationFiles);
    }
  }

  protected FilterConcatIntermediateFiles getNonMatedFilterConcatIntermediateFiles(MatchResult results, boolean paired, HashingRegion[] regions) throws IOException {
    Diagnostic.userLog("Processing " + results.size() + " hits for " + (paired ? "unmated " : "") + "reads");
    // for (MatchResult mr : results) {
    //   Diagnostic.log(mr.toString());
    // }

    if (mParams.useTopRandom()) {
      mSharedResources.getSingleEndTopRandom().initialize();
    }
    Diagnostic.progress("UnmatedInit: 1/1 Jobs Finished");

    final SimpleThreadPool stp = new SimpleThreadPool(mParams.numberThreads(), paired ? "UnmatedProcessing" : "Alignment", true);
    final MapQScoringReadBlocker blockerFirst = new MapQScoringReadBlockerSynch((int) mParams.buildFirstParams().reader().numberSequences(), mParams.outputParams().maxTopResults());

    final MapQScoringReadBlocker blockerSecond;
    if (paired) {
      blockerSecond = new MapQScoringReadBlockerSynch((int) mParams.buildFirstParams().reader().numberSequences(), mParams.outputParams().maxTopResults());
    } else {
      blockerSecond = null;
    }

    final String namePrefix = paired ? TEMP_SAM_UNMATED_ALIGNMENT_NAME : TEMP_SAM_ALIGNMENT_NAME;
    final File[] tempFiles = new File[regions.length];
    stp.enableBasicProgress(tempFiles.length);
    final ChunkPair[] chunks = findChunkBoundaries(regions, results);
    for (int i = 0; i < tempFiles.length; ++i) {
      tempFiles[i] = mParams.outputParams().resultStreamHandler().tempFile(namePrefix + i + FileUtils.GZ_SUFFIX);
      final OutputStream stream = FileUtils.createOutputStream(tempFiles[i], true);

      if (paired) {
        final PairedTempFileWriterImpl sw = new PairedTempFileWriterImpl(mParams,  mUnmappedTracker, mSharedResources);
        sw.initialiseUnmated(stream, blockerFirst, blockerSecond);
        stp.execute(new UnmatedWriterThread(sw, results, chunks[i].mChunkStart, chunks[i].mChunkEnd, regions[i], i));
      } else {
        final SingleEndTempFileWriter sw = new SingleEndTempFileWriter(mParams, mUnmappedTracker, mSharedResources);
        sw.initialiseAlignments(stream, blockerFirst);
        stp.execute(new AlignmentWriterThread(sw, results, chunks[i].mChunkStart, chunks[i].mChunkEnd, regions[i], i));
      }



    }
    stp.terminate();
    results.evacuateTheBuilding();

    //merge results
    if (paired) {
      Diagnostic.userLog("Merging unmated results");
    } else {
      Diagnostic.userLog("Merging alignment results");
    }
    final SingleEndTopRandomImplementation.HitRecord[] hitsToKeep;
    final NamesInterface templateNames;
    if (mParams.useTopRandom()) {
      hitsToKeep = mSharedResources.getSingleEndTopRandom().getRecords();
      templateNames = mSharedResources.names();
    } else {
      hitsToKeep = null;
      templateNames = null;
    }
    final File outFile;
    if (mParams.outputParams().bam()) {
      outFile = mParams.outputParams().resultStreamHandler().file(paired ? NgsOutputParams.UNMATED_BAM_FILE_NAME : NgsOutputParams.ALIGNMENTS_BAM_FILE_NAME);
    } else {
      outFile = mParams.outputParams().isCompressOutput()
              ? mParams.outputParams().resultStreamHandler().file((paired ? NgsOutputParams.UNMATED_SAM_FILE_NAME : NgsOutputParams.ALIGNMENTS_SAM_FILE_NAME) + FileUtils.GZ_SUFFIX)
              : mParams.outputParams().resultStreamHandler().file(paired ? NgsOutputParams.UNMATED_SAM_FILE_NAME : NgsOutputParams.ALIGNMENTS_SAM_FILE_NAME);
    }


    final FilterConcatIntermediateFiles unmatedIntFiles = filterConcatNonMated(blockerFirst, blockerSecond, tempFiles, hitsToKeep, templateNames, outFile);

    //SamSingleEndOutputProcessor.indexSamFile(mParams, outFile);

    if (mParams.useTopRandom()) {
      mSharedResources.getSingleEndTopRandom().finish();
    }
    return unmatedIntFiles;
  }

  /**
   * Finds chunk boundaries corresponding to <code>regions</code>. That is return the pair of indexes corresponding to the first match within the padded region and the first match after the region. If no results are within the padded region then the pair of indexes will correspond to the first result after the padded region, or the length of the results array if no such result exists.
   * @param regions regions to find chunk boundaries for
   * @param results results to allocate to regions
   * @return as described above.
   */
  static ChunkPair[] findChunkBoundaries(HashingRegion[] regions, MatchResult results) {
    final ChunkPair[] chunks = new ChunkPair[regions.length];
    final long[] chunkStart = new long[regions.length];
    final long[] chunkEnd = new long[regions.length];
    Arrays.fill(chunkStart, -1L);
    int regionIndexStart = 0;
    int regionIndexEnd = 0;
    for (long resultIndex = 0; resultIndex < results.size(); ++resultIndex) {
      int cmpStart = -1;
      //find first region containing result
      while (regionIndexStart < regions.length) {
        cmpStart = regions[regionIndexStart].isInPaddedRange(results.getTemplateId(resultIndex), results.getPosition(resultIndex));
        if (cmpStart <= 0) {
          //result is left or in region
          break;
        } else {
          //result is right of region, look for region containing result
          ++regionIndexStart;
        }
      }
      if (cmpStart < 0) {
        continue; //can never find region containing result. move on
      }
      if (regionIndexEnd < regionIndexStart) {
        // regions had no results in them, make their boundaries be zero length, from and to first result outside their region
        for (int i = regionIndexEnd; i < regionIndexStart; ++i) {
          chunkStart[i] = resultIndex;
          chunkEnd[i] = resultIndex;
        }
        regionIndexEnd = regionIndexStart;
      }
      //find first region past result
      while (regionIndexEnd < regions.length) {
        final int cmpEnd = regions[regionIndexEnd].isInPaddedRange(results.getTemplateId(resultIndex), results.getPosition(resultIndex));
        if (cmpEnd >= 0) {
          //result is in region (+1 case should not happen given above constraints)
          ++regionIndexEnd;
        } else {
          //result is left of region, first region past result
          break;
        }
      }

      for (int i = regionIndexStart; i < regionIndexEnd; ++i) {
        if (chunkStart[i] == -1) {
          chunkStart[i] = resultIndex;
        }
        if (chunkEnd[i] < resultIndex + 1) {
          chunkEnd[i] = resultIndex + 1; //chunkEnd is exclusive, and everything between regionIndexStart and regionIndexEnd includes current result
        }
      }
    }
    for (int i = regionIndexEnd; i < regions.length; ++i) {
      chunkStart[i] = results.size();
      chunkEnd[i] = results.size();
    }
    for (int i = 0; i < regions.length; ++i) {
      chunks[i] = new ChunkPair(chunkStart[i], chunkEnd[i]);
    }
    return chunks;
  }

  @Override
  public OutputProcessor threadClone(HashingRegion region) throws IOException {
    mRegions.add(region);
    return this;
  }

  protected void sortRegions() {
    Collections.sort(mRegions);
  }

  static class ChunkPair {
    final long mChunkStart; //0 based inclusive
    final long mChunkEnd; //0 based exclusive

    ChunkPair(long chunkStart, long chunkEnd) {
      mChunkStart = chunkStart;
      mChunkEnd = chunkEnd;
    }

    @Override
    public boolean equals(Object o) {
      if (this == o) {
        return true;
      }
      if (o == null || getClass() != o.getClass()) {
        return false;
      }
      final ChunkPair chunkPair = (ChunkPair) o;
      return mChunkEnd == chunkPair.mChunkEnd && mChunkStart == chunkPair.mChunkStart;
    }

    @Override
    public int hashCode() {
      final int result = (int) (mChunkStart ^ (mChunkStart >>> 32));
      return 31 * result + (int) (mChunkEnd ^ (mChunkEnd >>> 32));
    }

    @Override
    public String toString() {
      return "ChunkPair{"
              + "mChunkStart=" + mChunkStart
              + ", mChunkEnd=" + mChunkEnd
              + '}';
    }
  }

  private static final class AlignmentWriterThread extends AbstractAlignmentWriterThread {

    AlignmentWriterThread(AbstractTempFileWriter writer, MatchResult results, long chunkStart, long chunkEnd, HashingRegion region, int threadNumber) {
      super(writer, results, chunkStart, chunkEnd, region, threadNumber);
    }

    @Override
    protected void handleResult(int templateId, int encodedReadId, int position, boolean reverse) throws IOException {
      ((SingleEndTempFileWriter) mSamWriter).alignmentResult(encodedReadId, reverse, position);
    }
  }

  private static class UnmatedWriterThread extends AbstractAlignmentWriterThread {

    UnmatedWriterThread(AbstractTempFileWriter writer, MatchResult results, long chunkStart, long chunkEnd, HashingRegion region, int threadNumber) {
      super(writer, results, chunkStart, chunkEnd, region, threadNumber);
    }

    @Override
    protected void handleResult(int templateId, int encodedReadId, int position, boolean reverse) throws IOException {
      final int readId = PAIRED_END.decode(encodedReadId);
      final boolean first = PAIRED_END.isFirst(encodedReadId);
      ((PairedTempFileWriterImpl) mSamWriter).unmatedResult(readId, first , reverse, position);
    }

  }


}
