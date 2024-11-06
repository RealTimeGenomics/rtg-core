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
package com.rtg.position.output;

/**
 */
public final class OutputUtils {

  private OutputUtils() { } //prevent instantiation

  static class ScoreResult {
    private double mScore;
    private int mIdCount;

    void set(final double score, final int idCount) {
      mScore = score;
      mIdCount = idCount;
    }

    double score() {
      return mScore;
    }

    int idCount() {
      return mIdCount;
    }
  }

  /**
   * Compute best possible scoring based on BLOSUM62 matrix for a given gap.
   * @param res where the results are returned.
   * @param aGap length of one side of the gap.
   * @param bGap length of side of the gap.
   * @param wordSize number of residues in a word.
   * @param stepSize number of residues in a step.
   */
  static void gapScore(final ScoreResult res, final int aGap, final int bGap, final int wordSize, final int stepSize) {
    if (aGap > bGap) {
      gapScore(res, bGap, aGap, wordSize, stepSize);
      return;
    }
    assert aGap >= 0;
    assert bGap >= 0;
    final int subs;
    if (aGap == bGap) {
      assert aGap >= stepSize;
      if (aGap == stepSize) {
        subs = 1;
      } else {
        final int t = aGap - 2 * stepSize;
        final int w = t / wordSize;
        subs = w + 2;
      }
      final int idC = aGap - subs;
      res.set(idC * Blosum62.HIT - subs /* * Blosum62.MISS */, idC);
      //System.err.println(" buildLength=" + buildLength() + " scoreIncrement=" + mScoreIncrement + " score=" + mScore + " sc=" + sc + " GappedScoreTracks");
      return;
    }
    if (aGap <= stepSize) {
      //the deletion takes care of a single substitution
      subs = 0;
    } else {
      final int t = aGap - 2 * stepSize;
      final int w = t / wordSize;
      subs = w + 1;
    }
    final int indels = bGap - aGap;
    final int idC2 = aGap - subs;
    final double sc = idC2 * Blosum62.HIT - subs /* * Blosum62.MISS */ + Blosum62.GAP - indels/* * Blosum62.EXTEND */;
    //System.err.println(" buildLength=" + buildLength() + " scoreIncrement=" + mScoreIncrement + " score=" + mScore + " sc=" + sc + " GappedScoreTracks");
    res.set(sc, idC2);
  }
}
