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

package com.rtg.variant.bayes;

/**
 * @param <D> description type
 */
public class MockModel<D extends Description> extends Model<D> {
  /**
   * @param hypotheses description of the hypotheses.
   * @param statistics collects statistics about the evidence.
   * @param post dummy posteriors.
   */
  public MockModel(Hypotheses<D> hypotheses, final Statistics<?> statistics, final double[] post) {
    super(hypotheses, statistics, new NoAlleleBalance());
    if (post != null) {
      for (int i = 0; i < size(); i++) {
        mPosteriors[i] = arithmetic().prob2Poss(post[i]);
      }
    }
  }

}
