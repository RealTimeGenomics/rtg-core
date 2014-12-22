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
package com.rtg.ngs.tempstage;

import java.io.IOException;
import java.io.OutputStream;
import java.util.Comparator;

import com.rtg.alignment.ActionsHelper;
import com.rtg.alignment.AlignmentResult;
import com.rtg.ngs.NgsParams;
import com.rtg.ngs.ReadStatusTracker;
import com.rtg.ngs.SharedResources;
import com.rtg.ngs.blocking.MapQScoringReadBlocker;
import com.rtg.pairedend.ReadStatusListener;
import com.rtg.reader.ReadHelper;
import com.rtg.reader.SequencesReader;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * Write single end alignments to temp files
 *
 */
public class SingleEndTempFileWriter extends AbstractTempFileWriter {

  protected static final int MAPPED = ReadStatusTracker.UNMATED_FIRST | ReadStatusTracker.UNMATED_SECOND;

  private SmartTempFileWriter mBinarizableRecordWriter = null;
  private MapQScoringReadBlocker mUnmatedBlocker = null;

  private int mMaxScoreFailed = 0;
  private int mMaxScorePassed = 0;

  /**
   * Construct a new writer for the given read lane and template, writing
   * results to specified output directory.
   * @param param the parameters
   * @param listener a listener that will receive notifications of mated pairs
   * @param sharedResources resources shared between multiple instances of this class
   */
  public SingleEndTempFileWriter(NgsParams param, ReadStatusListener listener, SharedResources sharedResources) {
    super(listener, sharedResources, param);
  }

  @Override
  public void close() throws IOException {
    super.close();
    closeAlignments();
  }


  /**
   * Closes alignment output
   * @throws IOException if an IO exception occurs
   */
  public void closeAlignments() throws IOException {
    if (mBinarizableRecordWriter != null) {
      try {
        final int over = mBinarizableRecordWriter.getMaxCapacityUsed();
        //System.err.println("Overflow is " + over);
        //System.err.println("Dups: " +  mBinarizableRecordWriter.getDuplicateCount());
        Diagnostic.userLog("Alignment score filter (" + mMatedMaxMismatches + "): " + mMaxScorePassed + " passed, " + (mMaxScorePassed + mMaxScoreFailed) + " total.");

        Diagnostic.userLog("Reordering buffer used capacity of " + over + " records");
        Diagnostic.userLog("Duplicates detected during SAM writing: " + mBinarizableRecordWriter.getDuplicateCount());
      } finally {
        mBinarizableRecordWriter.close();
        mBinarizableRecordWriter = null;
      }
    }
  }


  /**
   * Tells the alignment writer where alignment results should be written.
   * Must be called before any calls to <code>alignmentResult</code>
   * @param alignmentsOut output stream to send unmated results to.
   * @param alignmentsBlocker blocker to keep track of best score per read
   * @exception IOException if an I/O error occurs
   */
  public void initialiseAlignments(final OutputStream alignmentsOut, final MapQScoringReadBlocker alignmentsBlocker) throws IOException {
    mUnmatedBlocker = alignmentsBlocker;
    mBinarizableRecordWriter = createSmartWriter(alignmentsOut);
    mTemplateId = -1;
  }

  /**
   * Align and write result
   * @param readId read identifier
   * @param rc true if match was on reverse strand
   * @param start start position of match
   * @return true if result was written, false if it was filtered
   * @throws IOException if an io exception occurs
   */
  public boolean alignmentResult(int readId, boolean rc, int start) throws IOException {
    final SequencesReader reader = mFirstReader;
    final byte[] read = ReadHelper.getRead(reader, readId);
    final int size = reader.length(readId);
    final int maxMismatches = mMatedMaxMismatches.getValue(size);
    final int maxScore = maxMismatches * mSubstitutionPenalty;
    final int[] matchActions = calculateEditDistance(read, read.length, start, rc, mMatedMaxMismatches, true, readId);
    mListener.addStatus(readId, ReadStatusTracker.UNMATED_COMPUTE_ALIGNMENT_FIRST);


    /* Double checker
    boolean dump = false;
    if (matchActions[ActionsHelper.ALIGNMENT_SCORE_INDEX] > maxScore) { // I.e. we termated early
      //final int newMaxScore = mMaxScore.getValue(size);
      final int newMaxScore = 1000000;
      matchActions = calculateEditDistance(read, read.length, start, rc, newMaxScore);
      final int newAS = matchActions[ActionsHelper.ALIGNMENT_SCORE_INDEX];
      if (newAS < maxScore) {
        Diagnostic.developerLog("Liar, liar, pants on fire. Alignment at " + start + " terminated early with maxscore of " + maxScore
                                + " but then found alignment " + (maxScore - newAS) + " better with score "
                                + newAS + " by shifting " + (matchActions[ActionsHelper.TEMPLATE_START_INDEX] + mTemplateOffset - start));
        dump=true;
      }
    }
    */

    if (!mClipRegion.isInRange(mTemplateId, matchActions[ActionsHelper.TEMPLATE_START_INDEX] + mTemplateOffset)) {
      // The alignment is out of range. Bail out.
      return true;
    }

    if (matchActions[ActionsHelper.ALIGNMENT_SCORE_INDEX] <= maxScore) {
      mMaxScorePassed++;
      mListener.addStatus(readId, ReadStatusTracker.UNMATED_ALIGN_SCORE_FIRST);
      final MapQScoringReadBlocker blocker = mUnmatedBlocker;
      if (!blocker.isBlocked2(readId, matchActions[ActionsHelper.ALIGNMENT_SCORE_INDEX])) {
        final AlignmentResult matchResult = new AlignmentResult(read, matchActions, mTemplate);
        matchResult.setIdentifyingInfo(true, rc);
        matchResult.setRemainingOutput(readId, (int) mTemplateId); //mSharedResources.names().name(mTemplateId));

        // Double checker: if (dump) Diagnostic.developerLog(matchResult.toString(mTemplateOffset));
        final BinaryTempFileRecord sam = matchResult.toRecord(false, null, mTemplateOffset, false, mLegacy);
        final boolean notDup = mBinarizableRecordWriter.addAlignmentHandleDuplicates(sam);
        if (notDup) {
          final int newScore = mUnmatedBlocker.increment(readId, matchResult.getScore());
          if (mSharedResources.getSingleEndTopRandom() != null && newScore != -1) {
            mSharedResources.getSingleEndTopRandom().update(readId, (int) mTemplateId, mTemplateOffset + matchResult.getStart(), matchResult.isReverse(), matchResult.getScore(), newScore);
          }
        }
        return true;
      }
    } else {
      mMaxScoreFailed++;
    }
    return false;
  }


  /**
   * Align and write result without applying top equals filtering (alignment score threshold is still applied).
   * @param readId read identifier
   * @param rc true if match was on reverse strand
   * @param start start position of match
   * @return the SAM record to be written, or null if outside the clip region
   * @throws IOException if an io exception occurs
   */
  public BinaryTempFileRecord alignmentResultUnfiltered(int readId, boolean rc, int start) throws IOException {
    final SequencesReader reader = mFirstReader;
    final byte[] read = ReadHelper.getRead(reader, readId);
    final int size = reader.length(readId);
    final int[] matchActions = calculateEditDistance(read, read.length, start, rc, mMatedMaxMismatches, true, readId);
    mListener.addStatus(readId, ReadStatusTracker.UNMATED_COMPUTE_ALIGNMENT_FIRST);

    if (!mClipRegion.isInRange(mTemplateId, matchActions[ActionsHelper.TEMPLATE_START_INDEX] + mTemplateOffset)) {
      // The alignment is out of range (we'll get this alignment in another thread). Bail out.
      return null;
    }

    final int maxMismatches = mMatedMaxMismatches.getValue(size);
    final int maxScore = maxMismatches * mSubstitutionPenalty;
    if (matchActions[ActionsHelper.ALIGNMENT_SCORE_INDEX] <= maxScore) {
      mMaxScorePassed++;
      mListener.addStatus(readId, MAPPED);
      mListener.addStatus(readId, ReadStatusTracker.UNMATED_ALIGN_SCORE_FIRST);
      final AlignmentResult matchResult = new AlignmentResult(read, matchActions, mTemplate);
      matchResult.setIdentifyingInfo(true, rc);
      matchResult.setRemainingOutput(readId, (int) mTemplateId); //mSharedResources.names().name(mTemplateId));

      final BinaryTempFileRecord record = matchResult.toRecord(false, null, mTemplateOffset, true, mLegacy);
      mBinarizableRecordWriter.addAlignmentHandleDuplicates(record);
      return record;
    } else {
      mMaxScoreFailed++;
      return null;
    }
  }

  @Override
  Comparator<BinaryTempFileRecord> getRecordComparator() {
    return new TempFileRecordComparator();
  }
}
