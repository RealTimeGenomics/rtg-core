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

import static com.rtg.ngs.ReadStatusTracker.UNMATED_FIRST;
import static com.rtg.ngs.ReadStatusTracker.UNMATED_SECOND;

import java.io.IOException;
import java.io.OutputStream;
import java.util.Comparator;

import com.rtg.alignment.ActionsHelper;
import com.rtg.alignment.AlignmentResult;
import com.rtg.index.hash.ngs.ReadEncoder;
import com.rtg.mode.DnaUtils;
import com.rtg.ngs.NgsParams;
import com.rtg.ngs.ReadStatusTracker;
import com.rtg.ngs.SharedResources;
import com.rtg.ngs.blocking.MapQScoringReadBlocker;
import com.rtg.pairedend.MatedHitInfo;
import com.rtg.pairedend.ReadStatusListener;
import com.rtg.reader.ReadHelper;
import com.rtg.reader.SequencesReader;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * Write paired reads to temp files
 *
 */
public class PairedTempFileWriterImpl extends AbstractTempFileWriter implements PairedTempFileWriter {

  // count statuses per mapping - if output per output line
  private int mMatedCachedAlignment = 0;
  private int mMatedNonCachedAlignment = 0;
  private int mMatedMaxScoreFailed = 0;
  private int mMatedMaxScorePassed = 0;


  private SmartTempFileWriter mBinarizableRecordWriter = null;
  private SmartTempFileWriter mBinarizableRecordUnmatedWriter = null;
  private MapQScoringReadBlocker mUnmatedBlockerFirst;
  private MapQScoringReadBlocker mUnmatedBlockerSecond;

  /**
   * Construct a new writer for given left and right read lanes and template, writing
   * results to specified output directory.
   * @param param the params
   * @param listener a listener that will receive notifications of mated pairs
   * @param sharedResources resources shared between multiple instances of this class
   */
  public PairedTempFileWriterImpl(NgsParams param, ReadStatusListener listener, SharedResources sharedResources) {
    super(listener, sharedResources, param);
  }

  private int[] mActionsBufferLeft = null;
  private void bufferActionsLeft(int[] actions) {
    if ((mActionsBufferLeft == null) || (mActionsBufferLeft.length < actions.length)) {
      mActionsBufferLeft = new int[actions.length];
    }
    System.arraycopy(actions, 0, mActionsBufferLeft, 0, actions.length);
  }
  private int[] mActionsBufferRight = null;
  private void bufferActionsRight(int[] actions) {
    if ((mActionsBufferRight == null) || (mActionsBufferRight.length < actions.length)) {
      mActionsBufferRight = new int[actions.length];
    }
    System.arraycopy(actions, 0, mActionsBufferRight, 0, actions.length);
  }

  @Override
  public boolean pairResultLeft(final MatedHitInfo matedHitInfo) throws IOException {
    if (mTemplate == null) {
      throw new RuntimeException("Call nextTemplateId first");
    }
    assert matedHitInfo.getReadId() < Integer.MAX_VALUE;
    final int readId = matedHitInfo.getReadId();
    mListener.addStatus(readId, ReadStatusTracker.MATED);

    byte[] read1 = null;
    final int score1;
    if (matedHitInfo.getAlignmentScoreLeft() != -1) {
      //Diagnostic.developerLog("Using cached left alignment score for read " + readId);
      ++mMatedCachedAlignment;
      score1 = matedHitInfo.getAlignmentScoreLeft();
    } else {
      ++mMatedNonCachedAlignment;
      read1 = ReadHelper.getRead(!matedHitInfo.isFirstRight() ? mFirstReader : mSecondReader, matedHitInfo.getReadId());
      bufferActionsLeft(calculateEditDistance(read1, read1.length, matedHitInfo.getTemplateStartLeft(), matedHitInfo.isReverseComplementLeft(), mMatedMaxMismatches, !matedHitInfo.isFirstRight(), readId));
      score1 = ActionsHelper.alignmentScore(mActionsBufferLeft);
      matedHitInfo.setAlignmentScoreLeft(score1);
    }

    // If leftmost side score *really* sucks (either exceeds alignment
    // score or enough to disqualify combo score on its own) pass this
    // status so sliding window collector can boot other pairings
    // containing it
    final int matchReadLength = !matedHitInfo.isFirstRight() ? mFirstReader.length(readId) : mSecondReader.length(readId);
    if ((score1 > mMatedMaxMismatches.getValue(matchReadLength) * mSubstitutionPenalty) || mSharedResources.getBlocker().isBlocked2(readId, score1)) {
      matedHitInfo.setLeftVeryPoor(true);
      return false;
    }

    byte[] read2 = null;
    final int score2;
    if (matedHitInfo.getAlignmentScoreRight() != -1) {
      //Diagnostic.developerLog("Using cached right alignment score for read " + readId);
      ++mMatedCachedAlignment;
      score2 = matedHitInfo.getAlignmentScoreRight();
    } else {
      ++mMatedNonCachedAlignment;
      read2 = ReadHelper.getRead(matedHitInfo.isFirstRight() ? mFirstReader : mSecondReader, matedHitInfo.getReadId());
      bufferActionsRight(calculateEditDistance(read2, read2.length, matedHitInfo.getTemplateStartRight(), matedHitInfo.isReverseComplementRight(), mMatedMaxMismatches, matedHitInfo.isFirstRight(), readId));
      score2 = ActionsHelper.alignmentScore(mActionsBufferRight);
      matedHitInfo.setAlignmentScoreRight(score2);
    }

    if (score2 == Integer.MAX_VALUE) {
      return false;
    }

    if (checkPairScores(readId, !matedHitInfo.isFirstRight(), score1, score2)) {
      AlignmentResult matchResult = matedHitInfo.getAlignmentLeft();
      if (matchResult == null) { //
        if (read1 == null) { // Score from another potential mating was used, still need to recalc alignment since that wasn't cached
          //Diagnostic.developerLog("Recalculating left alignment for read " + readId + " at pos " + matedHitInfo.mTemplateStartLeft);
          read1 = ReadHelper.getRead(!matedHitInfo.isFirstRight() ? mFirstReader : mSecondReader, matedHitInfo.getReadId());
          bufferActionsLeft(calculateEditDistance(read1, read1.length, matedHitInfo.getTemplateStartLeft(), matedHitInfo.isReverseComplementLeft(), mMatedMaxMismatches, !matedHitInfo.isFirstRight(), readId));
          if (mActionsBufferLeft[ActionsHelper.ALIGNMENT_SCORE_INDEX] == Integer.MAX_VALUE) {
            return false;
          }
        }
        matchResult = new AlignmentResult(read1, mActionsBufferLeft, mTemplate);
        matchResult.setIdentifyingInfo(!matedHitInfo.isFirstRight(), matedHitInfo.isReverseComplementLeft());
      }
      AlignmentResult mateResult = matedHitInfo.getAlignmentRight();
      if (mateResult == null) {
        if (read2 == null) { // Score from another potential mating was used, still need to recalc alignment since that wasn't cached
          //Diagnostic.developerLog("Recalculating right alignment for read " + readId + " at pos " + matedHitInfo.mTemplateStartRight);
          read2 = ReadHelper.getRead(matedHitInfo.isFirstRight() ? mFirstReader : mSecondReader, matedHitInfo.getReadId());
          bufferActionsRight(calculateEditDistance(read2, read2.length, matedHitInfo.getTemplateStartRight(), matedHitInfo.isReverseComplementRight(), mMatedMaxMismatches, matedHitInfo.isFirstRight(), readId));
        }
        if (mActionsBufferRight[ActionsHelper.ALIGNMENT_SCORE_INDEX] == Integer.MAX_VALUE) {
          return false;
        }
        mateResult = new AlignmentResult(read2, mActionsBufferRight, mTemplate);
        mateResult.setIdentifyingInfo(matedHitInfo.isFirstRight(), matedHitInfo.isReverseComplementRight());
      }

      matedHitInfo.setAlignments(matchResult, mateResult);

      pairResult(readId, matchResult, mateResult);
      return true;
    }

    return false;
  }

  @Override
  public void pairResultRight(final MatedHitInfo matedHitInfo) throws IOException {
    if (mTemplate == null) {
      throw new RuntimeException("Call nextTemplateId first");
    }
    final int readId = matedHitInfo.getReadId();
    final AlignmentResult matchResult = matedHitInfo.getAlignmentRight();
    final AlignmentResult mateResult = matedHitInfo.getAlignmentLeft();
    mListener.addStatus(readId, ReadStatusTracker.MATED);
    if (checkPairScores(readId, matchResult.isFirst(), matchResult.getScore(), mateResult.getScore())) {
      pairResult(readId, matchResult, mateResult);
    }
  }

  boolean checkPairScores(int readId, boolean matchFirst, int matchScore, int mateScore) throws IOException {
    final int comboScore = matchScore + mateScore;
    // In extreme cases the addition might overflow, because this is addition the
    // following test will catch the overflow.
    if (comboScore < 0) {
      throw new ArithmeticException("Overflow");
    }
    if (!mSharedResources.getBlocker().isBlocked2(readId, comboScore)) {
      final int matchReadLength = matchFirst ? mFirstReader.length(readId) : mSecondReader.length(readId);
      final int mateReadLength = !matchFirst ? mFirstReader.length(readId) : mSecondReader.length(readId);
      if (matchScore <= mMatedMaxMismatches.getValue(matchReadLength) * mSubstitutionPenalty && mateScore <= mMatedMaxMismatches.getValue(mateReadLength) * mSubstitutionPenalty) {
        ++mMatedMaxScorePassed;
        return true;
      } else {
        ++mMatedMaxScoreFailed;
      }
    }
    return false;
  }

  @Override
  Comparator<BinaryTempFileRecord> getRecordComparator() {
    return new PairedTempFileRecordComparator();
  }

  private void pairResult(int readId, AlignmentResult matchResult, AlignmentResult mateResult) throws IOException {
    // Only actually output if the current side is in the clip region.
    // Between this and the blocker we may get unbalanced records in an individual file,
    // but the final filterconcat makes everything consistent.
    // NOTE: if sliding window mated hit overflow occurs within the padding area,
    // there is the potential for unbalanced mated records in the final output.
    if (mClipRegion.isInRange(mTemplateId, matchResult.getStart() + mTemplateOffset)) {
      mListener.addStatus(readId, ReadStatusTracker.MATED_ALIGN_SCORE);

      final int comboScore = matchResult.getScore() + mateResult.getScore();
      matchResult.setRemainingOutput(readId, (int) mTemplateId);
      mateResult.setRemainingOutput(readId, (int) mTemplateId);
      final BinaryTempFileRecord rec = matchResult.toRecord(true, mateResult, mTemplateOffset, false, mLegacy);
      rec.setComboScore(comboScore);
      final boolean notDup = mBinarizableRecordWriter.addAlignmentHandleDuplicates(rec);
      // Only increment blocker once per pair
      if (matchResult.isFirst() && notDup) {
        final int newCount = mSharedResources.getBlocker().increment(readId, comboScore);
        if (mSharedResources.getPairedEndTopRandom() != null && newCount != -1) {
          mSharedResources.getPairedEndTopRandom().update(readId, (int) mTemplateId, matchResult.getStart() + mTemplateOffset, matchResult.isReverse(), mateResult.getStart() + mTemplateOffset, mateResult.isReverse(), comboScore, newCount);
        }
      }
    }
  }


  @Override
  public void close() throws IOException {
    super.close();
    closeMated();
    closeUnmated();
  }


  @Override
  public void closeMated() throws IOException {
    if (mBinarizableRecordWriter != null) {
      try (SmartTempFileWriter o = mBinarizableRecordWriter) {
        mBinarizableRecordWriter = null;
        final int over = o.getMaxCapacityUsed();
        //System.err.println("Overflow is " + over);
        //System.err.println("Dups: " +  mBinarizableRecordWriter.getDuplicateCount());
        final int scoreTotal = mMatedMaxScorePassed + mMatedMaxScoreFailed;
        Diagnostic.userLog("Alignment score filter (mismatch threshold of " + mMatedMaxMismatches + "): " + mMatedMaxScorePassed + " passed, " + scoreTotal + " total, " + (100.0 * mMatedMaxScorePassed / scoreTotal) + " %");

        final int alTotal = mMatedNonCachedAlignment + mMatedCachedAlignment;
        Diagnostic.developerLog("Alignments computed: " + mMatedNonCachedAlignment + " cached:" + mMatedCachedAlignment + " total, " + (100.0 * mMatedCachedAlignment / alTotal) + " %");

        Diagnostic.userLog("Reordering buffer used capacity of " + over + " records");
        Diagnostic.userLog("Duplicates detected during SAM writing: " + o.getDuplicateCount());
      }
    }
  }

  @Override
  @SuppressWarnings("try")
  public void closeUnmated() throws IOException {
    if (mBinarizableRecordUnmatedWriter != null) {
      try (SmartTempFileWriter ignored = mBinarizableRecordUnmatedWriter) {
        mBinarizableRecordUnmatedWriter = null;
      }
    }
  }

  @Override
  public void initialiseMated(OutputStream matedOut) throws IOException {
    mBinarizableRecordWriter = createSmartWriter(matedOut);
  }

  @Override
  public void initialiseUnmated(final OutputStream unmatedOut, final MapQScoringReadBlocker unmatedBlockerFirst, final MapQScoringReadBlocker unmatedBlockerSecond) throws IOException {
    mUnmatedBlockerFirst = unmatedBlockerFirst;
    mUnmatedBlockerSecond = unmatedBlockerSecond;
    mBinarizableRecordUnmatedWriter = createSmartWriter(unmatedOut);
    mTemplateId = -1;
  }


  @Override
  public boolean unmatedResult(int readId, boolean first, boolean rc, int start) throws IOException {
    final SequencesReader reader = first ? mFirstReader : mSecondReader;
    final MapQScoringReadBlocker blocker = first ? mUnmatedBlockerFirst : mUnmatedBlockerSecond;
    final byte[] read = ReadHelper.getRead(reader, readId);
    final int size = reader.length(readId);
    if ((start - mTemplateOffset < 0) || (start - mTemplateOffset > mTemplate.length)) {
      Diagnostic.developerLog("PROBLEM: start out of clipped template range."
                              + "\ntemplateid:" + mTemplateId + " start:" + start + " offset:" + mTemplateOffset + " trimmedlength:" + mTemplate.length
                              + "\n>readid" + readId + "\n" + DnaUtils.bytesToSequenceIncCG(read)
                              + "\n>readid" + readId + "\n" +  DnaUtils.bytesToSequenceIncCG(ReadHelper.getRead(!first ? mFirstReader : mSecondReader, readId)) + "\n");
    }
    final int[] matchActions = calculateEditDistance(read, read.length, start, rc, mUnmatedMaxMismatches, first, readId);
    mListener.addStatus(readId, first ? ReadStatusTracker.UNMATED_COMPUTE_ALIGNMENT_FIRST : ReadStatusTracker.UNMATED_COMPUTE_ALIGNMENT_SECOND);
    // if (readId == 1273) {
    //   Diagnostic.developerLog("read: " + readId + " seq:" + Utils.bytesToSequenceIncCG(read) + " rc:" + rc + " start:" + start + " offset:" + mTemplateOffset + " actions:" + ActionsHelper.toString(matchActions) + " shiftedto:" + matchActions[ActionsHelper.TEMPLATE_START_INDEX] + " region:" + mClipRegion);
    //   Diagnostic.developerLog("Templateid:" + mTemplateId + " template:" + Utils.bytesToSequenceIncCG(mTemplate));
    // }
    if (!mClipRegion.isInRange(mTemplateId, matchActions[ActionsHelper.TEMPLATE_START_INDEX] + mTemplateOffset)) {
      // The alignment is out of range. Bail out.
      return true;
    }

    final int maxMismatches = mUnmatedMaxMismatches.getValue(size);
    final int maxScore = maxMismatches * mSubstitutionPenalty;
    if (matchActions[ActionsHelper.ALIGNMENT_SCORE_INDEX] <= maxScore) {
      mListener.addStatus(readId, first ? ReadStatusTracker.UNMATED_ALIGN_SCORE_FIRST : ReadStatusTracker.UNMATED_ALIGN_SCORE_SECOND);
      if (!blocker.isBlocked2(readId, matchActions[ActionsHelper.ALIGNMENT_SCORE_INDEX])) {
        final AlignmentResult matchResult = new AlignmentResult(read, matchActions, mTemplate);
        matchResult.setIdentifyingInfo(first, rc);
        matchResult.setRemainingOutput(readId, (int) mTemplateId); //mSharedResources.names().name(mTemplateId));

        final BinaryTempFileRecord record = matchResult.toRecord(true, null, mTemplateOffset, false, mLegacy);
        final boolean notDup = mBinarizableRecordUnmatedWriter.addAlignmentHandleDuplicates(record);
        if (notDup) {
          final int encodedReadId;
          final int newCount;
          newCount = blocker.increment(readId, matchResult.getScore());
          if (first) {
            encodedReadId = ReadEncoder.PAIRED_FIRST.encode(readId);
          } else {
            encodedReadId = ReadEncoder.PAIRED_SECOND.encode(readId);
          }
          if (mSharedResources.getSingleEndTopRandom() != null && newCount != -1) {
            mSharedResources.getSingleEndTopRandom().update(encodedReadId, (int) mTemplateId, matchResult.getStart() + mTemplateOffset, matchResult.isReverse(), matchResult.getScore(), newCount);
          }
        }
        return true;
      }
    }
    return false;
  }

  /**
   * Write an unmated result without any top equals filtering (but still applying alignment score threshold)
   * @param readId read identifier
   * @param first true if read is first in pair
   * @param rc true if match was on reverse strand
   * @param start start position of match
   * @return the SAM record to be written, or null if outside the clip region
   * @throws IOException if an IO exception occurs
   */
  public BinaryTempFileRecord unmatedResultUnfiltered(int readId, boolean first, boolean rc, int start) throws IOException {
    final SequencesReader reader = first ? mFirstReader : mSecondReader;
    final byte[] read = ReadHelper.getRead(reader, readId);
    final int size = reader.length(readId);
    final int[] matchActions = calculateEditDistance(read, read.length, start, rc, mMatedMaxMismatches, first, readId);
    mListener.addStatus(readId, first ? ReadStatusTracker.UNMATED_COMPUTE_ALIGNMENT_FIRST : ReadStatusTracker.UNMATED_COMPUTE_ALIGNMENT_SECOND);
    if (!mClipRegion.isInRange(mTemplateId, matchActions[ActionsHelper.TEMPLATE_START_INDEX] + mTemplateOffset)) {
      // The alignment is out of range (we'll get this alignment in another thread). Bail out.
      return null;
    }


    final int maxMismatches = mMatedMaxMismatches.getValue(size);
    final int maxScore = maxMismatches * mSubstitutionPenalty;
    if (matchActions[ActionsHelper.ALIGNMENT_SCORE_INDEX] <= maxScore) {
      mListener.addStatus(readId, first ? ReadStatusTracker.UNMATED_ALIGN_SCORE_FIRST : ReadStatusTracker.UNMATED_ALIGN_SCORE_SECOND);
      mListener.addStatus(readId, first ? UNMATED_FIRST : UNMATED_SECOND);
      final AlignmentResult matchResult = new AlignmentResult(read, matchActions, mTemplate);
      matchResult.setIdentifyingInfo(first, rc);
      matchResult.setRemainingOutput(readId, (int) mTemplateId); //mSharedResources.names().name(mTemplateId));

      final BinaryTempFileRecord record = matchResult.toRecord(true, null, mTemplateOffset, true, mLegacy);
      mBinarizableRecordUnmatedWriter.addAlignmentHandleDuplicates(record);
      return record;
    } else {
      return null;
    }
  }
}
