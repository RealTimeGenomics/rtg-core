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

import java.util.Locale;

import com.rtg.mode.DnaUtils;

import junit.framework.TestCase;

/**
 */
public class CorrectingMutatorTest extends TestCase {
  public void testAllMutations() {
    String original = "act";
    String[] mutations = new String[]{
        "cct"
        , "gct"
        , "tct"
        , "aat"
        , "agt"
        , "att"
        , "aca"
        , "acc"
        , "acg"
    };
    int start = 0;
    int end = 3;
    check(original, mutations, start, end);
  }

  private void check(String original, String[] mutations, int start, int end) {
    CorrectingMutator cm = new CorrectingMutator();
    int i = 0;
    for (CorrectingMutator.SequenceBases mutant : cm.getMutations(new CorrectingMutator.BaseRead(DnaUtils.encodeString(original)), start, end)) {
      assertEquals(mutations[i], mutant.toString().toLowerCase(Locale.getDefault()));
      i++;
    }
    assertEquals(mutations.length, i);
  }

  public void testEndMutations() {
    String original = "act";
    String[] mutations = new String[]{
          "aat"
        , "agt"
        , "att"
        , "aca"
        , "acc"
        , "acg"
    };
    int start = 1;
    int end = 3;
    check(original, mutations, start, end);
  }

  public void testStartMutations() {
    String original = "act";
    String[] mutations = new String[]{
        "cct"
        , "gct"
        , "tct"
        , "aat"
        , "agt"
        , "att"
    };
    int start = 0;
    int end = 2;
    check(original, mutations, start, end);
  }
}
