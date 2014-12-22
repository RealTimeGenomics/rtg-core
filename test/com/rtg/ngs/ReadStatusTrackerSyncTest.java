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
package com.rtg.ngs;

import junit.framework.TestCase;

/**
 */
public class ReadStatusTrackerSyncTest extends TestCase {

  /**
   */
  public ReadStatusTrackerSyncTest(String name) {
    super(name);
  }

  public void testTracker() {
    final MapStatistics stats = new PairedEndMapStatistics(false, null);
    final ReadStatusTrackerSync trackerSync = new ReadStatusTrackerSync(7, stats);

    trackerSync.addStatus(0, ReadStatusTracker.MATCHED_SECOND);
    trackerSync.addStatus(0, ReadStatusTracker.MATCHED_FIRST);
    trackerSync.addStatus(0, ReadStatusTracker.MATED);

    assertEquals('d', trackerSync.getXCAttribute(0, false));

    trackerSync.addStatus(1, ReadStatusTracker.MATCHED_SECOND);
    trackerSync.addStatus(1, ReadStatusTracker.MATCHED_FIRST);
    trackerSync.addStatus(1, ReadStatusTracker.MATED);
    trackerSync.addStatus(1, ReadStatusTracker.MATED_ALIGN_SCORE);

    assertEquals('e', trackerSync.getXCAttribute(1, false));

    trackerSync.addStatus(2, ReadStatusTracker.MATCHED_FIRST);
    trackerSync.addStatus(2, ReadStatusTracker.MATCHED_SECOND);
    trackerSync.addStatus(2, ReadStatusTracker.UNMATED_ALIGN_SCORE_FIRST);

    assertEquals('C', trackerSync.getXCAttribute(2, true));

    trackerSync.addStatus(3, ReadStatusTracker.MATCHED_FIRST);
    trackerSync.addStatus(3, ReadStatusTracker.MATCHED_SECOND);
    trackerSync.addStatus(3, ReadStatusTracker.UNMATED_ALIGN_SCORE_FIRST);
    trackerSync.addStatus(3, ReadStatusTracker.UNMATED_COMPUTE_ALIGNMENT_FIRST);

    assertEquals('E', trackerSync.getXCAttribute(3, true));

    trackerSync.addStatus(4, ReadStatusTracker.MATCHED_FIRST);
    trackerSync.addStatus(4, ReadStatusTracker.MATCHED_SECOND);
    trackerSync.addStatus(4, ReadStatusTracker.UNMATED_ALIGN_SCORE_SECOND);
    trackerSync.addStatus(4, ReadStatusTracker.UNMATED_COMPUTE_ALIGNMENT_SECOND);

    assertEquals('E', trackerSync.getXCAttribute(4, false));

    trackerSync.addStatus(5, ReadStatusTracker.MATCHED_FIRST);
    trackerSync.addStatus(5, ReadStatusTracker.MATCHED_SECOND);
    trackerSync.addStatus(5, ReadStatusTracker.UNMATED_COMPUTE_ALIGNMENT_FIRST);

    assertEquals('D', trackerSync.getXCAttribute(5, true));

    trackerSync.addStatus(6, ReadStatusTracker.MATCHED_FIRST);
    trackerSync.addStatus(6, ReadStatusTracker.MATCHED_SECOND);
    trackerSync.addStatus(6, ReadStatusTracker.UNMATED_COMPUTE_ALIGNMENT_SECOND);

    assertEquals('D', trackerSync.getXCAttribute(6, false));

    trackerSync.addStatus(6, ReadStatusTracker.BLOCKED_SECOND);
    assertEquals('B', trackerSync.getXCAttribute(6, false));
    assertEquals('C', trackerSync.getXCAttribute(6, true));
    trackerSync.addStatus(6, ReadStatusTracker.BLOCKED_FIRST);
    assertEquals('B', trackerSync.getXCAttribute(6, true));
  }
}
