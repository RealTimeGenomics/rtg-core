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
