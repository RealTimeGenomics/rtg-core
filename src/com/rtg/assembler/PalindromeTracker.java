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

import java.util.HashSet;
import java.util.Set;

import com.rtg.assembler.graph.Graph;

/**
 */
public class PalindromeTracker {
  Set<Long> mSet = new HashSet<>();
  PalindromeTracker(Graph graph) {
    for (long i = 1; i < graph.numberContigs() + 1; i++) {
      if (graph.contigDeleted(i)) {
        continue;
      }
      if (GraphAlignment.isPalindrome(i, graph)) {
        mSet.add(i);
      }
    }
  }

  boolean isPalindrome(long contigId) {
    return mSet.contains(Math.abs(contigId));
  }
}
