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

import com.rtg.ngs.ReadStatusTracker.UnmappedStatus;
import com.rtg.util.TestUtils;

import junit.framework.TestCase;

/**
 */
public class ReadStatusTrackerTest extends TestCase {

  /**
   */
  public ReadStatusTrackerTest(String name) {
    super(name);
  }

  public void testStatusToString() {
    assertEquals("MATED_FIRST ", ReadStatusTracker.statusToString(ReadStatusTracker.MATED_FIRST));
    assertEquals("MATED_SECOND ", ReadStatusTracker.statusToString(ReadStatusTracker.MATED_SECOND));
    assertEquals("UNMATED_FIRST ", ReadStatusTracker.statusToString(ReadStatusTracker.UNMATED_FIRST));
    assertEquals("UNMATED_SECOND ", ReadStatusTracker.statusToString(ReadStatusTracker.UNMATED_SECOND));
    assertEquals("UNMAPPED_FIRST ", ReadStatusTracker.statusToString(ReadStatusTracker.UNMAPPED_FIRST));
    assertEquals("UNMAPPED_SECOND ", ReadStatusTracker.statusToString(ReadStatusTracker.UNMAPPED_SECOND));

    assertEquals("BLOCKED_FIRST ", ReadStatusTracker.statusToString(ReadStatusTracker.BLOCKED_FIRST));
    assertEquals("BLOCKED_SECOND ", ReadStatusTracker.statusToString(ReadStatusTracker.BLOCKED_SECOND));

    assertEquals("MATCHED_FIRST ", ReadStatusTracker.statusToString(ReadStatusTracker.MATCHED_FIRST));
    assertEquals("MATCHED_SECOND ", ReadStatusTracker.statusToString(ReadStatusTracker.MATCHED_SECOND));

    assertEquals("MATED_ALIGN_SCORE ", ReadStatusTracker.statusToString(ReadStatusTracker.MATED_ALIGN_SCORE));
    assertEquals("MATED ", ReadStatusTracker.statusToString(ReadStatusTracker.MATED));

    assertEquals("UNIQUELY_MAPPED_FIRST ", ReadStatusTracker.statusToString(ReadStatusTracker.UNIQUELY_MAPPED_FIRST));
    assertEquals("UNIQUELY_MAPPED_SECOND ", ReadStatusTracker.statusToString(ReadStatusTracker.UNIQUELY_MAPPED_SECOND));

    assertEquals("UNMATED_ALIGN_SCORE_FIRST ", ReadStatusTracker.statusToString(ReadStatusTracker.UNMATED_ALIGN_SCORE_FIRST));
    assertEquals("UNMATED_ALIGN_SCORE_SECOND ", ReadStatusTracker.statusToString(ReadStatusTracker.UNMATED_ALIGN_SCORE_SECOND));

    assertEquals("UNMATED_COMPUTE_ALIGNMENT_FIRST ", ReadStatusTracker.statusToString(ReadStatusTracker.UNMATED_COMPUTE_ALIGNMENT_FIRST));
    assertEquals("UNMATED_COMPUTE_ALIGNMENT_SECOND ", ReadStatusTracker.statusToString(ReadStatusTracker.UNMATED_COMPUTE_ALIGNMENT_SECOND));
  }

  public void testTracker() {
    final MapStatistics stats = new PairedEndMapStatistics(false, null);
    final ReadStatusTracker tracker = new ReadStatusTracker(9, stats);

    tracker.addStatus(0, ReadStatusTracker.MATCHED_SECOND);
    tracker.addStatus(0, ReadStatusTracker.MATCHED_FIRST);
    tracker.addStatus(0, ReadStatusTracker.MATED);

    assertEquals('d', tracker.getXCAttribute(0, false));

    tracker.addStatus(1, ReadStatusTracker.MATCHED_SECOND);
    tracker.addStatus(1, ReadStatusTracker.MATCHED_FIRST);
    tracker.addStatus(1, ReadStatusTracker.MATED);
    tracker.addStatus(1, ReadStatusTracker.MATED_ALIGN_SCORE);

    assertEquals('e', tracker.getXCAttribute(1, false));

    tracker.addStatus(2, ReadStatusTracker.MATCHED_FIRST);
    tracker.addStatus(2, ReadStatusTracker.MATCHED_SECOND);
    tracker.addStatus(2, ReadStatusTracker.UNMATED_ALIGN_SCORE_FIRST);

    assertEquals('C', tracker.getXCAttribute(2, true));

    tracker.addStatus(3, ReadStatusTracker.MATCHED_FIRST);
    tracker.addStatus(3, ReadStatusTracker.MATCHED_SECOND);
    tracker.addStatus(3, ReadStatusTracker.UNMATED_ALIGN_SCORE_FIRST);
    tracker.addStatus(3, ReadStatusTracker.UNMATED_COMPUTE_ALIGNMENT_FIRST);

    assertEquals('E', tracker.getXCAttribute(3, true));

    tracker.addStatus(4, ReadStatusTracker.MATCHED_FIRST);
    tracker.addStatus(4, ReadStatusTracker.MATCHED_SECOND);
    tracker.addStatus(4, ReadStatusTracker.UNMATED_ALIGN_SCORE_SECOND);
    tracker.addStatus(4, ReadStatusTracker.UNMATED_COMPUTE_ALIGNMENT_SECOND);

    assertEquals('E', tracker.getXCAttribute(4, false));

    tracker.addStatus(5, ReadStatusTracker.MATCHED_FIRST);
    tracker.addStatus(5, ReadStatusTracker.MATCHED_SECOND);
    tracker.addStatus(5, ReadStatusTracker.UNMATED_COMPUTE_ALIGNMENT_FIRST);

    assertEquals('D', tracker.getXCAttribute(5, true));

    tracker.addStatus(6, ReadStatusTracker.MATCHED_FIRST);
    tracker.addStatus(6, ReadStatusTracker.MATCHED_SECOND);
    tracker.addStatus(6, ReadStatusTracker.UNMATED_COMPUTE_ALIGNMENT_SECOND);

    assertEquals('D', tracker.getXCAttribute(6, false));

    tracker.addStatus(6, ReadStatusTracker.BLOCKED_SECOND);
    assertEquals('B', tracker.getXCAttribute(6, false));
    assertEquals('C', tracker.getXCAttribute(6, true));
    tracker.addStatus(6, ReadStatusTracker.BLOCKED_FIRST);
    assertEquals('B', tracker.getXCAttribute(6, true));

    assertEquals('A', tracker.getXCAttribute(7, true));
    assertEquals('A', tracker.getXCAttribute(7, false));
    tracker.addStatus(7, ReadStatusTracker.BLOCKED_FIRST);
    tracker.addStatus(7, ReadStatusTracker.BLOCKED_SECOND);
    assertEquals('B', tracker.getXCAttribute(7, true));
    assertEquals('B', tracker.getXCAttribute(7, false));

    tracker.addStatus(8, ReadStatusTracker.MATED);
    assertEquals('d', tracker.getXCAttribute(8, true));
    assertEquals('d', tracker.getXCAttribute(8, false));
    tracker.addStatus(8, ReadStatusTracker.MATED_ALIGN_SCORE);
    assertEquals('e', tracker.getXCAttribute(8, true));
    assertEquals('e', tracker.getXCAttribute(8, false));


  }

  public void testGetStatus() {
    final MapStatistics stats = new PairedEndMapStatistics(false, null);
    final ReadStatusTracker tracker = new ReadStatusTracker(7, stats);

    tracker.addStatus(0, ReadStatusTracker.MATCHED_SECOND);
    tracker.addStatus(0, ReadStatusTracker.MATCHED_FIRST);
    tracker.addStatus(0, ReadStatusTracker.MATED);

    assertTrue(tracker.getStatus(0, ReadStatusTracker.MATCHED_SECOND));
    assertTrue(tracker.getStatus(0, ReadStatusTracker.MATCHED_FIRST));
    assertTrue(tracker.getStatus(0, ReadStatusTracker.MATED));

    assertFalse(tracker.getStatus(0, ReadStatusTracker.BLOCKED_FIRST));
    assertFalse(tracker.getStatus(0, ReadStatusTracker.BLOCKED_SECOND));
    assertFalse(tracker.getStatus(0, ReadStatusTracker.MATED_FIRST));
    assertFalse(tracker.getStatus(0, ReadStatusTracker.MATED_SECOND));
    assertFalse(tracker.getStatus(0, ReadStatusTracker.UNMATED_COMPUTE_ALIGNMENT_FIRST));
    assertFalse(tracker.getStatus(0, ReadStatusTracker.UNMATED_COMPUTE_ALIGNMENT_SECOND));
    assertFalse(tracker.getStatus(0, ReadStatusTracker.UNMATED_ALIGN_SCORE_FIRST));
    assertFalse(tracker.getStatus(0, ReadStatusTracker.UNMATED_ALIGN_SCORE_SECOND));
    assertFalse(tracker.getStatus(0, ReadStatusTracker.UNMATED_FIRST));
    assertFalse(tracker.getStatus(0, ReadStatusTracker.UNMATED_SECOND));
    assertFalse(tracker.getStatus(0, ReadStatusTracker.UNMAPPED_FIRST));
    assertFalse(tracker.getStatus(0, ReadStatusTracker.UNMAPPED_SECOND));
    assertFalse(tracker.getStatus(0, ReadStatusTracker.UNIQUELY_MAPPED_FIRST));
    assertFalse(tracker.getStatus(0, ReadStatusTracker.UNIQUELY_MAPPED_SECOND));
  }

  public void testUnmappedStatusEnum() {
    TestUtils.testEnum(UnmappedStatus.class, "[LEFT_UNMAPPED, RIGHT_UNMAPPED, BOTH_UNMAPPED, SINGLE_END_UNMAPPED, MAPPED]");
  }
}
