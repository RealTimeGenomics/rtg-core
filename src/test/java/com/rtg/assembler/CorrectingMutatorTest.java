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

package com.rtg.assembler;

import java.util.Locale;

import com.rtg.mode.DnaUtils;

import junit.framework.TestCase;

/**
 */
public class CorrectingMutatorTest extends TestCase {
  public void testAllMutations() {
    final String original = "act";
    final String[] mutations = new String[]{
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
    final int start = 0;
    final int end = 3;
    check(original, mutations, start, end);
  }

  private void check(String original, String[] mutations, int start, int end) {
    final CorrectingMutator cm = new CorrectingMutator();
    int i = 0;
    for (CorrectingMutator.SequenceBases mutant : cm.getMutations(new CorrectingMutator.BaseRead(DnaUtils.encodeString(original)), start, end)) {
      assertEquals(mutations[i], mutant.toString().toLowerCase(Locale.getDefault()));
      ++i;
    }
    assertEquals(mutations.length, i);
  }

  public void testEndMutations() {
    final String original = "act";
    final String[] mutations = new String[]{
          "aat"
        , "agt"
        , "att"
        , "aca"
        , "acc"
        , "acg"
    };
    final int start = 1;
    final int end = 3;
    check(original, mutations, start, end);
  }

  public void testStartMutations() {
    final String original = "act";
    final String[] mutations = new String[]{
        "cct"
        , "gct"
        , "tct"
        , "aat"
        , "agt"
        , "att"
    };
    final int start = 0;
    final int end = 2;
    check(original, mutations, start, end);
  }
}
