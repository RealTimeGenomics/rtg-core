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
public class InvertedEnvironment extends EnvironmentImplementation {

  /**
   * @param env environment to be inverted.
   */
  public InvertedEnvironment(final EnvironmentImplementation env) {
    super(env.maxShift(), env.mTemplate, env.absoluteTemplatePosition(0), env.mRead, env.mQuality);
  }

  /**
   * @param env environment to be inverted.
   * @param absReadPos the expected absolute template position of the start of the read (in the original non-inverted environment).
   */
  public InvertedEnvironment(final EnvironmentImplementation env, final int absReadPos) {
    super(env.maxShift(), env.mTemplate, absReadPos, env.mRead, env.mQuality);
  }

  @Override
  public double quality(int index) {
    return super.quality(readLength() - index - 1) ;
  }

  @Override
  public byte template(int index) {
    return super.template(readLength() - index - 1);
  }

  @Override
  public byte read(int index) {
    return super.read(readLength() - index - 1) ;
  }

  @Override
  public int absoluteTemplatePosition(int index) {
    return super.absoluteTemplatePosition(readLength() - index - 1);
  }
}
