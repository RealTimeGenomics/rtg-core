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
