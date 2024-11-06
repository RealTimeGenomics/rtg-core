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
package com.rtg.pairedend;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Properties;

import com.rtg.ngs.NgsParams;
import com.rtg.ngs.SharedResources;
import com.rtg.ngs.tempstage.PairedTempFileWriter;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.intervals.ReferenceRegions;
import com.rtg.util.machine.MachineOrientation;

/**
 * Receive matching reports from the two halves of a read and generate
 * pairings of those reads.
 *
 * Task 110. Implementation of sliding window.  Needs to ensure pairs
 * are put out as two SAM records in the correct template
 * order. Probably can ignore change in <code>templateId</code> except as signal to
 * flush and to pass to PairedTempFileWriter.
 *
 */
public class SlidingWindowCollector extends AbstractSlidingWindowCollector<HitInfo> {
  private final PairedTempFileWriter mAlignmentWriter; // thing to write

  private final ArrayList<MatedHitInfo>[] mMatedReadsWindow; // the sliding window, contains hits whose mate has already been written
  private final int[] mMatedReadsWindowInUse;

  // statistics counters
  private long mMatingCount = 0;
  private long mMaxMatedHitsExceededCount = 0;

  /**
   * Creates a new <code>SlidingWindowCollector</code> with a
   * window computed from left and right read length and
   * <code>alignmentWriter</code> to send results to.
   *
   * @param maxFragmentLength maximum distance between ends of read pair
   * @param minFragmentLength minimum distance between ends of read pair
   * @param pairOrientation required mating orientation
   * @param alignmentWriter a <code>PairedAlignmentWriter</code>
   * @param sharedResources resources shared between multiple instances of this class
   * @param calibrationRegions regions which reads target, or null if whole genome is covered
   */
  public SlidingWindowCollector(int maxFragmentLength, int minFragmentLength, MachineOrientation pairOrientation, PairedTempFileWriter alignmentWriter, SharedResources sharedResources, ReferenceRegions calibrationRegions) {
    super(maxFragmentLength, minFragmentLength, pairOrientation, sharedResources, calibrationRegions);

    if (alignmentWriter == null) {
      throw new NullPointerException();
    }

    @SuppressWarnings("unchecked")
    final ArrayList<MatedHitInfo>[] matedReadsWindow = (ArrayList<MatedHitInfo>[]) new ArrayList<?>[mWindowSize];
    mMatedReadsWindow = matedReadsWindow;
    mMatedReadsWindowInUse = new int[mWindowSize];
    for (int i = 0; i < mWindowSize; ++i) {
      mMatedReadsWindow[i] = new ArrayList<>();
    }

    mAlignmentWriter = alignmentWriter;
  }

  /**
   * Is used for tests to build shared resources from parameters.
   *
   * @param maxFragmentLength maximum distance between ends of read pair
   * @param minFragmentLength minimum distance between ends of read pair
   * @param pairOrientation required mating orientation
   * @param alignmentWriter a <code>PairedAlignmentWriter</code>
   * @param params to build shared resources
   * @throws IOException if shared resources cannot be built
   *
   */
  SlidingWindowCollector(int maxFragmentLength, int minFragmentLength, MachineOrientation pairOrientation, PairedTempFileWriter alignmentWriter, NgsParams params) throws IOException {
    this(maxFragmentLength, minFragmentLength, pairOrientation, alignmentWriter, SharedResources.generateSharedResources(params), params.outputParams().calibrateRegions());
  }

  private MatedHitInfo getMatedHitInfo(int i) {
    final MatedHitInfo ret;
    if (mMatedReadsWindowInUse[i] == getMaxHitsPerPosition() - 1) {
      ++mMaxMatedHitsExceededCount;
      if (mMaxMatedHitsExceededCount < 5) {
        Diagnostic.userLog("Max mated hits per position exceeded at template: " + mReferenceId + " templateStart: " + (mReadsWindow[i].size() > 0 ? "" + mReadsWindow[i].get(0).mTemplateStart : "unknown"));
      }
    }
    if (mMatedReadsWindowInUse[i] == getMaxHitsPerPosition()) {
      ret = null; //sorry, no vacancy.
    } else {
      if (mMatedReadsWindowInUse[i] < mMatedReadsWindow[i].size()) {
        ret = mMatedReadsWindow[i].get(mMatedReadsWindowInUse[i]);
      } else {
        ret = new MatedHitInfo();
        mMatedReadsWindow[i].add(ret);
      }
      mMatedReadsWindowInUse[i]++;
    }
    return ret;
  }

  private void returnMatedHitInfo(int i) {
    mMatedReadsWindowInUse[i]--;
  }

  @Override
  boolean checkPair(HitInfo hit, HitInfo mate) throws IOException {
    // also check for < window size for +64 case
    // check start positions within window size - for extra 64 case...
    final int mateSlot = windowPosition(mate.mTemplateStart);
    final MatedHitInfo matedhit = getMatedHitInfo(mateSlot);

    if (matedhit != null) {
      matedhit.setValues(hit.readId(), mate.first(), hit.reverseComplement(), hit.templateStart(), mate.reverseComplement(), mate.templateStart()
      );
      matedhit.setAlignments(hit.alignment(), mate.alignment());

      final boolean pairOk = mAlignmentWriter.pairResultLeft(matedhit);

      /* Cache alignment scores and potentially alignment results for each hit */
      if (matedhit.mAlignmentScoreLeft != -1) {
        hit.setAlignment(matedhit.mAlignmentScoreLeft, matedhit.mAlignmentLeft);
      }
      if (matedhit.mAlignmentScoreRight != -1) {
        mate.setAlignment(matedhit.mAlignmentScoreRight, matedhit.mAlignmentRight);
      }

      if (!pairOk) {
        final boolean leftSucks = matedhit.mLeftIsVeryPoor;
        returnMatedHitInfo(mateSlot); // Return the MatedHitInfo to the pool
        if (leftSucks) {
          return false;
        }
      }
      ++mMatingCount;
    }
    return true;
  }

  private void outputExistingMates(final ArrayList<MatedHitInfo> hits, int size) throws IOException {
    for (int i = 0; i < size; ++i) {
      final MatedHitInfo hit = hits.get(i);
      // prior mapping - just output
//      System.err.println("Right Mating: " + hit + " <==> ");
      mAlignmentWriter.pairResultRight(hit);
      //mAlignmentWriter.pairResult(hit.mReadId, hit.mFirstRight, hit.mReverseComplementRight, hit.mTemplateStartRight, hit.mReverseComplementLeft, hit.mTemplateStartLeft);
    }
  }

  @Override
  void flushToPosition(final int newStart) throws IOException {
    assert integrity();
//    System.err.println("F2P: " + mCurrentReferencePosition + " : " + newStart);
    while (mCurrentReferencePosition < newStart) {
      // set current to new start
      ++mCurrentReferencePosition;

      final int windowIndex = windowPosition(mCurrentReferencePosition);
//      System.err.println("flushing position " + mCurrentReferencePosition);
      final ArrayList<HitInfo> hits = mReadsWindow[windowIndex];
      if (mReadsWindowInUse[windowIndex] > 0) {

        // output left side of mates (this can find mates at the same position, so do this before outputting right side of mates)
        findNewMates(hits, mReadsWindowInUse[windowIndex]);

        // clear processed hits
        clearHits(hits, mReadsWindowInUse[windowIndex]);
        mReadsWindowInUse[windowIndex] = 0;
      } else {
        mReadsWindowInUse[windowIndex] = 0;
      }

      // output right side of mates
      if (mMatedReadsWindowInUse[windowIndex] > 0) {
        outputExistingMates(mMatedReadsWindow[windowIndex], mMatedReadsWindowInUse[windowIndex]);
        mMatedReadsWindowInUse[windowIndex] = 0;
      }
    }
  }


  @Override
  void writerNextTemplateId(long templateId) throws IOException {
    mAlignmentWriter.nextTemplateId(templateId);
  }

  @Override
  public Properties getStatistics() {
    final Properties stats = super.getStatistics();
    stats.setProperty("matings", Long.toString(mMatingCount));
    stats.setProperty("max_mated_hits_exceeded", Long.toString(mMaxMatedHitsExceededCount));
    return stats;
  }

  @Override
  HitInfo createHitInfo() {
    return new HitInfo();
  }
}
