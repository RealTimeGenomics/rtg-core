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
package com.rtg.ngs.tempstage;

import static com.rtg.ngs.ReadStatusTracker.UNMATED_FIRST;
import static com.rtg.ngs.ReadStatusTracker.UNMATED_SECOND;

import java.io.IOException;
import java.io.OutputStream;
import java.util.Comparator;

import com.rtg.alignment.ActionsHelper;
import com.rtg.alignment.AlignmentResult;
import com.rtg.ngs.NgsParams;
import com.rtg.ngs.ReadStatusTracker;
import com.rtg.ngs.SharedResources;
import com.rtg.ngs.blocking.ReadBlocker;
import com.rtg.pairedend.ReadStatusListener;
import com.rtg.pairedend.UnfilteredHitInfo;
import com.rtg.reader.ReadHelper;
import com.rtg.util.IntegerOrPercentage;


/**
 * Writes unfiltered records to temp files
 */
public class UnfilteredTempFileWriter extends AbstractTempFileWriter {

  private SmartTempFileWriter mBinarizableRecordUnmatedWriter = null;

  protected final ReadBlocker mFreqBlockerLeft;
  protected final ReadBlocker mFreqBlockerRight;

  /**
   * Construct a new writer.
   * @param listener a listener that will receive notifications of mated pairs
   * @param sharedResources resources shared between multiple instances of this class
   * @param params {@link NgsParams} for current run
   * @param left left read frequency blocker
   * @param right right read frequency blocker
   */
  public UnfilteredTempFileWriter(ReadStatusListener listener, SharedResources sharedResources, NgsParams params, ReadBlocker left, ReadBlocker right) {
    super(listener, sharedResources, params);
    mFreqBlockerLeft = left;
    mFreqBlockerRight = right;
  }

  private int[] mActionsBufferHit = null;
  private void bufferActionsHit(int[] actions) {
    if ((mActionsBufferHit == null) || (mActionsBufferHit.length < actions.length)) {
      mActionsBufferHit = new int[actions.length];
    }
    System.arraycopy(actions, 0, mActionsBufferHit, 0, actions.length);
  }
  private int[] mActionsBufferMate = null;
  private void bufferActionsMate(int[] actions) {
    if ((mActionsBufferMate == null) || (mActionsBufferMate.length < actions.length)) {
      mActionsBufferMate = new int[actions.length];
    }
    System.arraycopy(actions, 0, mActionsBufferMate, 0, actions.length);
  }

  @Override
  Comparator<BinaryTempFileRecord> getRecordComparator() {
    return new PairedTempFileRecordComparator();
  }

  /**
   * Initialise the unmated output
   * @param unmatedOut stream to output unmated to
   */
  public void initialiseUnmated(final OutputStream unmatedOut) {
    mBinarizableRecordUnmatedWriter = createSmartWriter(unmatedOut);
    mTemplateId = -1;
  }

  /**
   * Checks whether an unmated hit passes alignment score thresholds
   * @param hit an {@link UnfilteredHitInfo} containing the hit coordinates
   * @return true if result score is acceptable
   * @throws IOException if an IO error occurs
   */
  public boolean checkUnmatedScore(UnfilteredHitInfo hit) throws IOException {
    final int unmatedMaxScoreHit = getMaxScoreValue(mUnmatedMaxMismatches, hit.readId(), hit.first());
    if (hit.score() != -1) {
      return hit.score() <= unmatedMaxScoreHit;
    }

    //mMatedNonCachedAlignment++;
    final byte[] readHit = ReadHelper.getRead(hit.first() ? mFirstReader : mSecondReader, hit.readId());
    bufferActionsHit(calculateEditDistance(readHit, readHit.length, hit.templateStart(), hit.reverseComplement(), mUnmatedMaxMismatches, hit.first(), hit.readId()));
    mListener.addStatus(hit.readId(), hit.first() ? ReadStatusTracker.UNMATED_COMPUTE_ALIGNMENT_FIRST : ReadStatusTracker.UNMATED_COMPUTE_ALIGNMENT_SECOND);

    final int scoreHit = ActionsHelper.alignmentScore(mActionsBufferHit);
    if (scoreHit > unmatedMaxScoreHit) {
      hit.setAlignment(scoreHit, null);
      return false;
    }

    hit.setUnmatedOk(true);
    hit.setAlignment(scoreHit, new AlignmentResult(readHit, mActionsBufferHit, mTemplate));
    return true;
  }

  private int getMaxScoreValue(IntegerOrPercentage maxIntOrPercent, long readId, boolean first) throws IOException {
    return maxIntOrPercent.getValue(first ? mFirstReader.length(readId) : mSecondReader.length(readId)) * mSubstitutionPenalty;
  }

  /**
   * Checks whether an mated hit pair passes alignment score thresholds
   * @param hit an {@link UnfilteredHitInfo} containing the hit coordinates
   * @param mate another {@link UnfilteredHitInfo} containing the hit coordinates for the mate
   * @return true if the pair was considered a good pairing
   * @throws IOException if an IO error occurs
   */
  public boolean checkScores(UnfilteredHitInfo hit, UnfilteredHitInfo mate) throws IOException {
    final int scoreHit;
    final int scoreMate;
    byte[] readHit = null;
    byte[] readMate = null;
    final IntegerOrPercentage maxMaxMismatches = mMatedMaxMismatches.compareTo(mUnmatedMaxMismatches) < 0 ? mUnmatedMaxMismatches : mMatedMaxMismatches;
    boolean updateHitAlignment = false;
    boolean updateMateAlignment = false;
    if (hit.score() != -1) {
      //Diagnostic.developerLog("Using cached left alignment score for read " + readId);
      //mMatedCachedAlignment++;
      scoreHit = hit.score();
    } else {
      //mMatedNonCachedAlignment++;
      readHit = ReadHelper.getRead(hit.first() ? mFirstReader : mSecondReader, hit.readId());
      bufferActionsHit(calculateEditDistance(readHit, readHit.length, hit.templateStart(), hit.reverseComplement(), maxMaxMismatches, hit.first(), hit.readId()));
      mListener.addStatus(hit.readId(), hit.first() ? ReadStatusTracker.UNMATED_COMPUTE_ALIGNMENT_FIRST : ReadStatusTracker.UNMATED_COMPUTE_ALIGNMENT_SECOND);

      scoreHit = ActionsHelper.alignmentScore(mActionsBufferHit);
      updateHitAlignment = true;
    }

    //stop here if hit score above BOTH thresholds?
    final int matedMaxScoreHit = getMaxScoreValue(mMatedMaxMismatches, hit.readId(), hit.first());
    final int unmatedMaxScoreHit = getMaxScoreValue(mUnmatedMaxMismatches, hit.readId(), hit.first());
    if (scoreHit > matedMaxScoreHit && scoreHit > unmatedMaxScoreHit) {
      if (updateHitAlignment) {
        hit.setAlignment(scoreHit, null);
      }
      return false;
    }

    if (scoreHit <= matedMaxScoreHit) {
      if (mate.score() != -1) {
        //Diagnostic.developerLog("Using cached left alignment score for read " + readId);
        //mMatedCachedAlignment++;
        scoreMate = mate.score();
      } else {
        //mMatedNonCachedAlignment++;
        readMate = ReadHelper.getRead(mate.first() ? mFirstReader : mSecondReader, hit.readId());
        bufferActionsMate(calculateEditDistance(readMate, readMate.length, mate.templateStart(), mate.reverseComplement(), maxMaxMismatches, mate.first(), mate.readId()));
        mListener.addStatus(hit.readId(), mate.first() ? ReadStatusTracker.UNMATED_COMPUTE_ALIGNMENT_FIRST : ReadStatusTracker.UNMATED_COMPUTE_ALIGNMENT_SECOND);

        scoreMate = ActionsHelper.alignmentScore(mActionsBufferMate);
        updateMateAlignment = true;
      }

      final int matedMaxScoreMate = getMaxScoreValue(mMatedMaxMismatches, mate.readId(), mate.first());
      final int unmatedMaxScoreMate = getMaxScoreValue(mUnmatedMaxMismatches, mate.readId(), mate.first());
      if (scoreMate <= matedMaxScoreMate) {
        hit.setMatedOk(true);
        mate.setMatedOk(true);
        if (updateHitAlignment) {
          hit.setAlignment(scoreHit, new AlignmentResult(readHit, mActionsBufferHit, mTemplate));
        }
        if (updateMateAlignment) {
          mate.setAlignment(scoreMate, new AlignmentResult(readMate, mActionsBufferMate, mTemplate));
        }
        return true;
      }

      if (scoreMate <= unmatedMaxScoreMate) {
        mate.setUnmatedOk(true);
        if (updateMateAlignment) {
          mate.setAlignment(scoreMate, new AlignmentResult(readMate, mActionsBufferMate, mTemplate));
        }
      } else {
        if (updateMateAlignment) {
          mate.setAlignment(scoreMate, null);
        }
      }
    }

    if (scoreHit <= unmatedMaxScoreHit) {
      hit.setUnmatedOk(true);
      if (updateHitAlignment) {
        hit.setAlignment(scoreHit, new AlignmentResult(readHit, mActionsBufferHit, mTemplate));
      }
    } else {
      if (updateHitAlignment) {
        hit.setAlignment(scoreHit, null);
      }
      return false;
    }
    return true;
  }

  /**
   * Write an unmated result without any top equals filtering
   * Only <code>matedOk</code> or <code>unmatedOk</code> hits should be passed.
   * @param hit the hit to write out as unmated
   * @throws IOException if an IO exception occurs
   */
  public void unmatedResultUnfiltered(UnfilteredHitInfo hit) throws IOException {
    assert hit.getMatedOk() || hit.getUnmatedOk();
    assert hit.alignment() != null;

    if (!mClipRegion.isInRange(mTemplateId, hit.alignment().getStart() + mTemplateOffset)) {
      // The alignment is out of range (we'll get this alignment in another thread). Bail out.
      return;
    }

    hit.alignment().setIdentifyingInfo(hit.first(), hit.reverseComplement());
    hit.alignment().setRemainingOutput(hit.readId(), (int) mTemplateId); //mSharedResources.names().name(mTemplateId));

    final BinaryTempFileRecord record = hit.alignment().toRecord(true, null, mTemplateOffset, true, mLegacy);
    if (hit.getMatedOk()) {
      record.setUnfilteredMated(true);
    }

    final ReadBlocker rb = hit.first() ? mFreqBlockerLeft : mFreqBlockerRight;
    synchronized (rb) {
      if (rb.isBlocked(hit.readId())) {
        return;
      }
      mListener.addStatus(hit.readId(), hit.first() ? ReadStatusTracker.UNMATED_ALIGN_SCORE_FIRST : ReadStatusTracker.UNMATED_ALIGN_SCORE_SECOND);
      mListener.addStatus(hit.readId(), hit.first() ? UNMATED_FIRST : UNMATED_SECOND);
      mBinarizableRecordUnmatedWriter.addAlignmentHandleDuplicates(record);
      rb.increment(hit.readId()); //because this hit already passed the alignment score thresholds, it's ok to increment rb here
    }
  }

  @Override
  public void close() throws IOException {
    super.close();
    closeUnmated();
  }

  void closeUnmated() throws IOException {
    if (mBinarizableRecordUnmatedWriter != null) {
      mBinarizableRecordUnmatedWriter.close();
      mBinarizableRecordUnmatedWriter = null;
    }
  }
}
