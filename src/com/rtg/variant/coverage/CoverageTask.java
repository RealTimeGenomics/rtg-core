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
package com.rtg.variant.coverage;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.rtg.launcher.ParamsTask;
import com.rtg.mode.DnaUtils;
import com.rtg.reader.ReaderUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.sam.CircularBufferMultifileSinglePassReaderWindow;
import com.rtg.sam.SamReadingContext;
import com.rtg.sam.SamUtils;
import com.rtg.sam.ThreadedMultifileIteratorWrapper;
import com.rtg.tabix.TabixIndexer;
import com.rtg.tabix.UnindexableDataException;
import com.rtg.util.MathUtils;
import com.rtg.util.Populator;
import com.rtg.util.SingletonPopulatorFactory;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ParallelProgress;
import com.rtg.util.diagnostic.Timer;
import com.rtg.util.intervals.RangeList;
import com.rtg.util.intervals.RangeList.RangeData;
import com.rtg.util.intervals.ReferenceRanges;
import com.rtg.variant.bayes.multisample.ChunkInfo;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;

/**
 * The new and improved coverage module.
 */
public class CoverageTask extends ParamsTask<CoverageParams, CoverageStatistics> {
  // Scale chosen so that IH up to 11 can be exactly represented, but
  // still small enough that only ~15 bits are used for the fraction.
  // Testing shows this is actually more accurate than IEEE double
  // arithmetic for the situation of interest
  private static final double SCALE = 16.0 * 9.0 * 5 * 7 * 11;
  private static final double INV_SCALE = 1.0 / SCALE;
  private long[] mChunkCovPrev;
  private long[] mChunkCov;
  private int[] mIH1Prev;
  private int[] mIH1;
  private int[] mIHgt1Prev;
  private int[] mIHgt1;
  private CircularBufferMultifileSinglePassReaderWindow<CoverageReaderRecord> mCircularBuffer;
  private ChunkInfo mInfo;
  private int mCurrentReferenceLength;
  private Long mReferenceSequenceIndex;
  private byte[] mReferenceBytes;
  private byte[] mPrevReferenceBytes;
  private int mChunkStart;
  private int mChunkEnd;
  private int mChunkNumber;
  private ParallelProgress mPP = null;
  private Map<String, Long> mReferenceNames = null;

  private ThreadedMultifileIteratorWrapper<CoverageReaderRecord> mWrapper;

  /**
   * @param params parameters for coverage run
   * @param reportStream stream to send summary to
   * @param stats object to collect summary statistics
   */
  public CoverageTask(CoverageParams params, OutputStream reportStream, CoverageStatistics stats) {
    super(params, reportStream, stats, null);
  }


  String formatMetaForBed(final List<String> metaStrings) {
    final StringBuilder sb = new StringBuilder();
    int stringsOutput = 0;
    for (final String metaString : metaStrings) {
      if (metaString == null || "".equals(metaString.trim())) {
        continue;
      }
      if (stringsOutput > 0) {
        sb.append(",");
      }
      sb.append(metaString);
      stringsOutput++;
    }
    return sb.toString();
  }

  @Override
  protected void exec() throws IOException {
    final SamRecordCounter recCounts = new SamRecordCounter();
    final OutputStream os = mParams.bedStream();
    final CoverageWriter coverageWriter = new CoverageWriter(os, mParams);
    coverageWriter.init();
    final SAMFileHeader uberHeader = SamUtils.getUberHeader(mParams.mapped(), mParams.ignoreIncompatibleSamHeaders(), null);

    final SequencesReader reference = mParams.genome() == null ? null : mParams.genome().reader();
    if (reference != null) {
      SamUtils.checkUberHeaderAgainstReference(reference, uberHeader, mParams.ignoreIncompatibleSamHeaders());
      mReferenceNames = ReaderUtils.getSequenceNameMap(reference);
    } else {
      Diagnostic.warning("No reference supplied - unable to determine regions of unknown nucleotides.");
    }

    //This list will maintain the sort order of stats output.
    //TODO: if we want multithreading by running multiple refs at once, this name list should be broken into per reference lists.
    final Set<String> outputRegionNames = new LinkedHashSet<>();
    try {
      final SingletonPopulatorFactory<CoverageReaderRecord> pf = new SingletonPopulatorFactory<>(new CoverageReaderRecordPopulator());
      final SamReadingContext context = new SamReadingContext(mParams.mapped(), mParams.ioThreads(), mParams.filterParams(), uberHeader);
      final ReferenceRanges<String> ranges = context.referenceRanges();
      mWrapper = new ThreadedMultifileIteratorWrapper<>(context, pf);
      for (final SAMSequenceRecord r : uberHeader.getSequenceDictionary().getSequences()) {
        if (r.getSequenceLength() > 0) {
          final RangeList<String> rs = ranges.get(r.getSequenceName());
          if (rs != null) {
            mWrapper.setSequenceId(r.getSequenceIndex());
            processReference(coverageWriter, r, recCounts, rs, outputRegionNames);
          }
        }
      }
    } finally {
      coverageWriter.close();
      if (mWrapper != null) {
        mWrapper.close();
      }
    }
    recCounts.reportCounts();
    if (mParams.blockCompressed() && mParams.outputIndex()) {
      final Timer indexing = new Timer("CoverageIndex");
      indexing.start();
      final File file = mParams.outFile();
      try {
        if (mParams.bedOutput() || mParams.bedgraphOutput()) {
          new TabixIndexer(file, new File(file.getParent(), file.getName() + TabixIndexer.TABIX_EXTENSION)).saveBedIndex();
        } else {
          new TabixIndexer(file, new File(file.getParent(), file.getName() + TabixIndexer.TABIX_EXTENSION)).saveTsvIndex();
        }
      } catch (final UnindexableDataException e) {
        Diagnostic.warning("Cannot produce TABIX index for: " + file + ": " + e.getMessage());
      }
      indexing.stop();
      indexing.log();
    }
    mStatistics.setOutputRegionNames(outputRegionNames);

    //this is really a statistics output message, but easier to detect here at the moment...
    if (reference == null) {
      Diagnostic.userLog("Reference genome not specified, meaning no non-N counts.");
    }
  }


  /**
   * @param coverageWriter writer for coverage values
   * @param r sam sequence record for information about reference
   * @param recCounts record counter
   * @param rangeList null to process entire reference, or a RangeList providing ranges of interest.
   * @param nameContainsSet a set of output region names seen so far
   * @throws IOException if an exception occurs while reading or writing
   */
  private void processReference(CoverageWriter coverageWriter, SAMSequenceRecord r, SamRecordCounter recCounts, RangeList<String> rangeList, Set<String> nameContainsSet) throws IOException {
    final List<RangeList.RangeData<String>> ranges = rangeList.getRangeList();
    if (ranges.size() == 0) { //no ranges were specified for this reference, so bail out
      return;
    }
    mPP = new ParallelProgress(r.getSequenceName());

    setReferenceSequence(r, ranges.get(0).getStart(), ranges.get(ranges.size() - 1).getEnd());

    int rangeIndex = 0;
    RangeData<String> range = ranges.get(rangeIndex);
    mStatistics.setRange(range);
    assert mInfo.start() == range.getStart();
    addOriginalRangeToSortedNameList(range, nameContainsSet);

    final int smooth = mParams.smoothing();
    try {
      int lastValStartPos = range.getStart();
      int lastVal = -1;
      boolean refHasCoverage = false;

      final CoverageSmoothingWindow csw = new CoverageSmoothingWindow(mInfo.start(), mInfo.end(), smooth);
      csw.initWindow();

      //main coverage loop.
      for (int currentTemplatePosition = mInfo.start(); currentTemplatePosition < mInfo.end(); currentTemplatePosition++) {
        if (currentTemplatePosition == range.getEnd()) { //reached the (exclusive) end of a range. do things necessary at the end of a range BEFORE processing the base at this position.
          final boolean noCoverage = lastVal == 0 && lastValStartPos == range.getStart();
          if (!mParams.tsvOutput() && !(mParams.onlyMappedRegions() && noCoverage)) { //since this is the end of a range, output a 'region of different coverage' line.
            coverageWriter.setBedLabel(formatMetaForBed(range.getMeta()));
            coverageWriter.finalCoverageRegion(r.getSequenceName(), lastValStartPos, currentTemplatePosition, lastVal);
          }

          //get a new range
          rangeIndex++;
          if (rangeIndex >= ranges.size()) {
            break;
          }
          range = ranges.get(rangeIndex);
          mStatistics.setRange(range);
          addOriginalRangeToSortedNameList(range, nameContainsSet);

          //deal with gaps and reset variables
          lastVal = -1;
          lastValStartPos = range.getStart();
        }

        //now deal with the base at this template position.

        if (currentTemplatePosition < range.getStart()) { //we're outside a range. keep track of the sliding window for consistency's sake.
          csw.step();
        } else { //we're within a range.
          final int currCov = csw.currentCoverage();

          if (mParams.tsvOutput()) { //tsv outputs something at every position within a range.
            coverageWriter.setBedLabel(formatMetaForBed(range.getMeta()));
            coverageWriter.finalCoveragePosition(r.getSequenceName(), currentTemplatePosition, getIH1ForPosition(currentTemplatePosition), getIHgt1ForPosition(currentTemplatePosition), getCoverageForPosition(currentTemplatePosition) * INV_SCALE);
          }

          if (lastVal != -1 && currCov != lastVal) { //the coverage value at this position is different from the previous
            if (!mParams.tsvOutput()) {
              coverageWriter.setBedLabel(formatMetaForBed(range.getMeta()));
              coverageWriter.finalCoverageRegion(r.getSequenceName(), lastValStartPos, currentTemplatePosition, lastVal);
            }
            lastValStartPos = currentTemplatePosition;
          }
          lastVal = currCov;
          csw.step();

          //update statistics for this base.
          final double nonSmoothCov = getCoverageForPosition(currentTemplatePosition) * INV_SCALE;
          final byte base = getBaseForPosition(currentTemplatePosition);
          mStatistics.updateCoverageHistogram(nonSmoothCov, mParams.genome() != null && base == DnaUtils.UNKNOWN_RESIDUE, mParams.minimumCoverageForBreadth());

          if (!refHasCoverage && currCov > 0) {
            refHasCoverage = true;
          }
        }
      }

      mStatistics.setRange(null); //flush the last range out of statistics

      recCounts.incrementCounts(mCircularBuffer);

      if ((!mParams.onlyMappedRegions() || refHasCoverage) && mCurrentReferenceLength > 0 && !mParams.tsvOutput()) { //output the last coverage value data.
        coverageWriter.setBedLabel(formatMetaForBed(range.getMeta()));
        coverageWriter.finalCoverageRegion(r.getSequenceName(), lastValStartPos, mInfo.end(), lastVal);
      }

      mPP.updateProgress(100);
      Diagnostic.progress("Finished: " + r.getSequenceName());

    } finally {
      mCircularBuffer.close();
    }
  }

  private void addOriginalRangeToSortedNameList(RangeData<String> range, Set<String> sortedRegionNames) {
    for (final RangeList.RangeData<String> or : range.getOriginalRanges()) {
      assert or.getMeta().size() == 1;
      sortedRegionNames.add(or.getMeta().get(0));
    }
  }

  /**
   * Keeps track of the coverage smoothing window
   */
  private final class CoverageSmoothingWindow {
    private int mNumInDaSum = 0;
    private long mCoverageSum = 0;
    private int mPos;
    private final int mEndPos;
    private final int mSmoothWindowSize;

    /**
     *
     * @param startPos start position, zero based inclusive
     * @param endPos end position, zero based exclusive
     * @param smoothWindowSize size of the window to smooth over
     */
    CoverageSmoothingWindow(int startPos, int endPos, int smoothWindowSize) {
      mPos = startPos;
      mEndPos = endPos;
      mSmoothWindowSize = smoothWindowSize;
    }

    private int currentCoverage() {
      return (int) MathUtils.round((mCoverageSum * INV_SCALE) / mNumInDaSum);
    }

    private void initWindow() throws IOException {
      for (int i = mPos; i < mPos + mSmoothWindowSize + 1 && i < mEndPos; i++) {
        mCoverageSum += getCoverageForPosition(i);
        mNumInDaSum++;
      }
    }

    private void step() throws IOException {
      if (mPos + mSmoothWindowSize + 1 < mInfo.end()) {
        mCoverageSum += getCoverageForPosition(mPos + mSmoothWindowSize + 1);
        mNumInDaSum++;
      }
      if (mPos - mSmoothWindowSize >= mInfo.start()) {
        mCoverageSum -= getCoverageForPosition(mPos - mSmoothWindowSize);
        mNumInDaSum--;
      }
      mPos++;
    }
  }

  private void setReferenceSequence(SAMSequenceRecord r, int restrictionStart, int restrictionEnd) {
    mCurrentReferenceLength = r.getSequenceLength();
    if (mParams.genome() != null) {
      mReferenceSequenceIndex = mReferenceNames.get(r.getSequenceName());
    }
    final int chunkSize = mParams.smoothing() < (mParams.chunkSize() / 2) ? mParams.chunkSize() : (mParams.smoothing() * 2 + 2);
    mInfo = new ChunkInfo(r.getSequenceLength(), r.getSequenceName(), chunkSize, restrictionStart, restrictionEnd, mParams.execThreads(), 1000);
    if (mCircularBuffer != null) {
      mCircularBuffer.close();
    }

    final CoverageReaderRecordPopulator populator = new CoverageReaderRecordPopulator();
    mCircularBuffer = new CircularBufferMultifileSinglePassReaderWindow<>(mWrapper, populator, r.getSequenceIndex(), mInfo.start(), Integer.MAX_VALUE);
    mChunkCov = null;
    mChunkStart = -1;
    mChunkEnd = -1;
    mChunkNumber = 0;
  }

  private byte getBaseForPosition(int sequencePosition) {
    if (mParams.genome() == null) {
      return DnaUtils.UNKNOWN_RESIDUE;
    }
    final int currChunkPos = sequencePosition - mChunkStart;
    if (currChunkPos < 0) {
      assert -currChunkPos < mInfo.chunkSize();
      return mPrevReferenceBytes[mInfo.chunkSize() + currChunkPos]; //currChunkPos is -ve
    } else {
      assert currChunkPos < mInfo.chunkSize();
      return mReferenceBytes[currChunkPos];
    }
  }

  private int getValueForPosition(int sequencePosition, int[] arr, int[] prevArr) {
    final int currChunkPos = sequencePosition - mChunkStart;
    if (currChunkPos < 0) {
      assert -currChunkPos < mInfo.chunkSize();
      return prevArr[mInfo.chunkSize() + currChunkPos]; //currChunkPos is -ve
    } else {
      assert currChunkPos < mInfo.chunkSize();
      return arr[currChunkPos];
    }
  }

  private int getIH1ForPosition(int sequencePosition) {
    return getValueForPosition(sequencePosition, mIH1, mIH1Prev);
  }

  private int getIHgt1ForPosition(int sequencePosition) {
    return getValueForPosition(sequencePosition, mIHgt1, mIHgt1Prev);
  }

  private long getCoverageForPosition(int sequencePosition) throws IOException {
    while (sequencePosition >= mChunkEnd) {
      mPP.updateProgress(mInfo.percent(mChunkEnd));
      if (!loadNextChunk()) {
        throw new ArrayIndexOutOfBoundsException(sequencePosition);
      }
    }
    final int currChunkPos = sequencePosition - mChunkStart;
    if (currChunkPos < 0) {
      assert -currChunkPos < mInfo.chunkSize();
      return mChunkCovPrev[mChunkCovPrev.length + currChunkPos]; //currChunkPos is -ve
    } else {
      assert currChunkPos < mInfo.chunkSize();
      return mChunkCov[currChunkPos];
    }
  }

  private boolean loadNextChunk() throws IOException {
    if (mChunkNumber < mInfo.numberChunks()) {
      mChunkCovPrev = mChunkCov;
      mChunkCov = new long[mInfo.chunkSize()];
      if (mParams.tsvOutput()) {
        mIH1Prev = mIH1;
        mIH1 = new int[mInfo.chunkSize()];
        mIHgt1Prev = mIHgt1;
        mIHgt1 = new int[mInfo.chunkSize()];
      }
      if (mChunkNumber > 0) {
        mCircularBuffer.flush(mChunkStart, mChunkEnd);
      }
      mChunkStart = mChunkNumber * mInfo.chunkSize() + mInfo.start();
      mChunkEnd = Math.min(mChunkStart + mInfo.chunkSize(), mInfo.end());
      mPrevReferenceBytes = mReferenceBytes;
      if (mReferenceSequenceIndex != null) {
        mReferenceBytes = new byte[mInfo.chunkSize()];
        mParams.genome().reader().read(mReferenceSequenceIndex, mReferenceBytes, mChunkStart, mChunkEnd - mChunkStart);
      }
      mChunkNumber++;
      final Iterator<CoverageReaderRecord> it = mCircularBuffer.recordsOverlap(mChunkStart, mChunkEnd);
      while (it.hasNext()) {
        final CoverageReaderRecord crr = it.next();
        for (int j = 0; j < crr.getCoverageBitSet().length(); j++) {
          if (crr.getCoverageBitSet().get(j)) {
            final int index = crr.getStart() + j - mChunkStart;
            if (index >= 0 && index < mChunkCov.length) {
              mChunkCov[index] += MathUtils.round(crr.getCoverageMultiplier() * SCALE);
              if (mParams.tsvOutput()) {
                if (crr.getIH() == 1) {
                  mIH1[index]++;
                } else {
                  mIHgt1[index]++;
                }
              }
            }
          }
        }
      }
      return true;
    }
    return false;
  }


  private static class CoverageReaderRecordPopulator implements Populator<CoverageReaderRecord> {
    @Override
    public CoverageReaderRecord overflow(int position, int length) {
      throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public CoverageReaderRecord populate(SAMRecord source) {
      return new CoverageReaderRecord(source, 0);
    }
  }

  private static class SamRecordCounter {
    private int mValidRecords = 0;
    private int mInvalidRecords = 0;
    private int mFilteredRecords = 0;
    private int mTossedRecords = 0;

    void incrementCounts(CircularBufferMultifileSinglePassReaderWindow<?> cbmrw) {
      mValidRecords += cbmrw.getValidRecordsCount();
      mInvalidRecords += cbmrw.getInvalidRecordsCount();
      mFilteredRecords += cbmrw.getFilteredRecordsCount();
      mTossedRecords += cbmrw.getTossedRecordCount();
    }

    void reportCounts() {
      final String invalidRecordsWarning = mInvalidRecords + " records skipped because of SAM format problems.";
      if (mInvalidRecords > 0) {
        Diagnostic.warning(invalidRecordsWarning);
      } else {
        Diagnostic.userLog(invalidRecordsWarning);
      }
      if (mFilteredRecords > 0) {
        Diagnostic.userLog(mFilteredRecords + " records skipped due to input filtering criteria");
      }
      if (mTossedRecords > 0) {
        Diagnostic.userLog(mTossedRecords + " records skipped in extreme coverage regions");
      }
      Diagnostic.userLog(mValidRecords + " records processed");
    }
  }
}
