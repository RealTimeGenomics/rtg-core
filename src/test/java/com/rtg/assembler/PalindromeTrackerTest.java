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

package com.rtg.assembler;

import com.rtg.assembler.graph.Graph;

import junit.framework.TestCase;

/**
 */
public class PalindromeTrackerTest extends TestCase {
  public void testPalindromeTracker() {
    final Graph g = GraphMapCliTest.makeGraph(2, new String[] {"ACGT", "AACT", "ATTGTTAACAAT", "ACGCGA", "ACCGGT"}, new long[][] {});
    final PalindromeTracker tracker = new PalindromeTracker(g);
    assertTrue(tracker.isPalindrome(1));
    assertTrue(tracker.isPalindrome(3));
    assertTrue(tracker.isPalindrome(5));
    assertFalse(tracker.isPalindrome(2));
    assertFalse(tracker.isPalindrome(4));

  }
}
