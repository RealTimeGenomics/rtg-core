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

package com.rtg.variant.realign;

/**
 * An environment with the read and the template reversed.
 * Note that <code>templatePosition(i)</code> moves DOWN the template
 * as <code>i</code> increase, so <code>templateStart()</code> actually
 * points to the highest position on the template that the read is
 * expected to align with.
 *
 * NOTE: it would be nice to do this inversion using EnvironmentDecorator.
 * However, that doesn't work for the <code>toString(sb)</code> method,
 * because it calls other methods in the same class - inheritance handles
 * this correctly, whereas delegation does not give the desired effect.
 *
 */
public class InvertCgTemplateEnvironment extends EnvironmentCombined {

  /**
   * @param env environment to be inverted.
   */
  public InvertCgTemplateEnvironment(final EnvironmentCombined env) {
    super(env.mSamEnv, env.mSamEnv.start(), env.maxShift(), env.mTemEnv);
  }

  @Override
  public byte template(int index) {
    return super.template(ScoreMatrixCG.EXPECTED_LENGTH - index - 1);
  }

  @Override
  public int absoluteTemplatePosition(int index) {
    return super.absoluteTemplatePosition(ScoreMatrixCG.EXPECTED_LENGTH - index - 1);
  }


}
