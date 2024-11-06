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
