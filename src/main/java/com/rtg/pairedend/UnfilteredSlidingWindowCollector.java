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

import com.rtg.ngs.SharedResources;
import com.rtg.ngs.tempstage.UnfilteredTempFileWriter;
import com.rtg.util.intervals.ReferenceRegions;
import com.rtg.util.machine.MachineOrientation;


/**
 * A Sliding Window Collector implementation for Unfiltered paired end
 */
public class UnfilteredSlidingWindowCollector extends AbstractSlidingWindowCollector<UnfilteredHitInfo> {

  private final UnfilteredTempFileWriter mWriter;

  /**
   * Constructor for thread cloning.
   * @param maxInsertSize maximum insert size
   * @param minInsertSize minimum insert size
   * @param pairOrientation required mating orientation
   * @param alignmentWriter alignment writer
   * @param sharedResources the shared resources between threads
   * @param calibrationRegions regions which reads target, or null if whole genome is covered
   */
  public UnfilteredSlidingWindowCollector(int maxInsertSize, int minInsertSize, MachineOrientation pairOrientation, UnfilteredTempFileWriter alignmentWriter, SharedResources sharedResources, ReferenceRegions calibrationRegions) {
    super(maxInsertSize, minInsertSize, pairOrientation, sharedResources, calibrationRegions);
    mWriter = alignmentWriter;
  }


  @Override
  UnfilteredHitInfo createHitInfo() {
    return new UnfilteredHitInfo();
  }

  @Override
  boolean checkPair(UnfilteredHitInfo hit, UnfilteredHitInfo mate) throws IOException {
    //check the scores, set the ok flagness
    return mWriter.checkScores(hit, mate);
  }

  @Override
  void flushToPosition(final int newStart) throws IOException {
    assert integrity();
    //System.err.println("F2P: " + mCurrentTemplatePosition + " : " + newStart);
    while (mCurrentReferencePosition < newStart) {
      // set current to new start
      ++mCurrentReferencePosition;

      final int windowIndex = windowPosition(mCurrentReferencePosition);

      //System.err.println("flushing position " + mCurrentTemplatePosition);
      if (mReadsWindowInUse[windowIndex] > 0) {
        final ArrayList<UnfilteredHitInfo> hits = mReadsWindow[windowIndex];

        findNewMates(hits, mReadsWindowInUse[windowIndex]);

        //output all hits at this location that are ok due to mating (previously set by findNewMates)
        // - also checks if it's ok even if unmated.
        for (int i = 0; i < mReadsWindowInUse[windowIndex]; ++i) {
          final UnfilteredHitInfo hit = hits.get(i);
          if (hit.score() == -1) {
            //ask alignment writer to calculate edit distance & see if this hit passes unmated threshold
            mWriter.checkUnmatedScore(hit);
          }
          if (hit.getMatedOk() || hit.getUnmatedOk()) {
            mWriter.unmatedResultUnfiltered(hit);
          }
        }

        // clear processed hits
        clearHits(hits, mReadsWindowInUse[windowIndex]);
      }
      mReadsWindowInUse[windowIndex] = 0;
    }
  }

  @Override
  void writerNextTemplateId(long templateId) throws IOException {
    mWriter.nextTemplateId(templateId);
  }
}
