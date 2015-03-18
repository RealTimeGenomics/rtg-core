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
package com.rtg.pairedend;

import java.io.IOException;
import java.util.ArrayList;

import com.rtg.ngs.SharedResources;
import com.rtg.ngs.tempstage.UnfilteredTempFileWriter;
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
   */
  public UnfilteredSlidingWindowCollector(int maxInsertSize, int minInsertSize, MachineOrientation pairOrientation, UnfilteredTempFileWriter alignmentWriter, SharedResources sharedResources) {
    super(maxInsertSize, minInsertSize, pairOrientation, sharedResources);
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
      mCurrentReferencePosition++;

      final int windowIndex = windowPosition(mCurrentReferencePosition);

      //System.err.println("flushing position " + mCurrentTemplatePosition);
      if (mReadsWindowInUse[windowIndex] > 0) {
        final ArrayList<UnfilteredHitInfo> hits = mReadsWindow[windowIndex];

        findNewMates(hits, mReadsWindowInUse[windowIndex]);

        //output all hits at this location that are ok due to mating (previously set by findNewMates)
        // - also checks if it's ok even if unmated.
        for (int i = 0; i < mReadsWindowInUse[windowIndex]; i++) {
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
