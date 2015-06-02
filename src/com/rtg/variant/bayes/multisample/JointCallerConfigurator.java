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

package com.rtg.variant.bayes.multisample;

import java.io.IOException;

import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantStatistics;

/**
 * @param <V> the type of statistics object used by this joitn caller.
 */
public interface JointCallerConfigurator<V extends VariantStatistics> {

  /**
   * Creates a configuration for the joint caller
   * @param params params with which to create the configuration
   * @param statistics overall statistics associated with the joint caller
   * @return an <code>AbstractJointCallerConfiguration</code>
   * @throws IOException in the event of an <code>IOException</code> occurring.
   */
  AbstractJointCallerConfiguration getConfig(final VariantParams params, V statistics) throws IOException;
}
