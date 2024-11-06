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
package com.rtg.protein;

import com.rtg.alignment.ActionsHelper;
import com.rtg.mode.ProteinScoringMatrix;

/**
 * Helper class for different scoring in {@link ProteinAlignmentResult}
 */
final class ScoringHelper {

  private ScoringHelper() { }

  private static final double INVERSE_LOG2 = 1.0 / Math.log(2.0);
  /**
   * Convert a SLIM alignment score (which is a negated raw score)
   * into a BLAST E-score.  We calculate the effective read length
   * using the new BLAST formula described in the
   * <a href="http://www.ch.embnet.org/software/blast_help.html">BLAST documentation</a>.
   *
   * Since our genome database usually contains very long sequences,
   * we use its actual size rather than calculating the effective size
   * of each sequence.
   *
   * @param alignScore the alignment score is the negated raw score
   * @param readLen the number of residues in the read
   * @param databaseLen the total number of residues in the database
   * @param matrix protein matrix to be used
   * @return the BLAST E-score
   */
  static double computeEScore(int alignScore, int readLen, double databaseLen, ProteinScoringMatrix matrix) {
    final double k = matrix.getK();
    final double h = matrix.getH();
    final double lambda = matrix.getLambda();
    // our alignScore is already negated...
    final double effectiveReadLen = Math.max(1.0, readLen + lambda * alignScore / h);
    return k * effectiveReadLen * databaseLen * Math.exp(lambda * alignScore);
  }

  /**
   * Convert a RTG alignment score (which is a negated raw score)
   * into a BLAST bit score.
   *
   * @param alignScore the alignment score is the negated raw score
   * @param matrix the matrix used to compute the score
   * @return the BLAST bit score
   */
  static double computeBitScore(int alignScore, ProteinScoringMatrix matrix) {
    final double lambda = matrix.getLambda();
    final double logK = matrix.getLogK();
    return (lambda * -alignScore - logK) * INVERSE_LOG2;
  }

  static int percentage(final int value, final int total) {
    if (total == 0) {
      return 0;
    }
    return (int) (100 * value / (double) total + 0.5);
  }

  static int percentId(final int[] actions) {
    return percentage(ActionsHelper.matchCount(actions), ActionsHelper.actionsCount(actions));
  }
}
