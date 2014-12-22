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

import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.Integrity;


/**
 * An environment that delegates to another environment.
 * This is a base class that leaves the environment unchanged.
 * It is intended to be used as the superclass of real environment decorators.
 *
 */
public class EnvironmentDecorator implements Environment, Integrity {

  /** the environment that we delegate to */
  final Environment mEnv;

  /**
   * @param env environment to be decorated.
   */
  EnvironmentDecorator(final Environment env) {
    mEnv = env;
  }

  @Override
  public int absoluteTemplatePosition(int index) {
    return mEnv.absoluteTemplatePosition(index);
  }

  @Override
  public int readLength() {
    return mEnv.readLength();
  }

  @Override
  public int maxShift() {
    return mEnv.maxShift();
  }

  @Override
  public double quality(int index) {
    return mEnv.quality(index);
  }

  @Override
  public byte read(int index) {
    return mEnv.read(index);
  }

  @Override
  public byte template(int index) {
    return mEnv.template(index);
  }

  @Override
  public int templateLength() {
    return mEnv.templateLength();
  }

  @Override
  public String toString() {
    return mEnv.toString();
  }

  @Override
  public boolean globalIntegrity() {
    return Exam.globalIntegrity(mEnv);
  }

  @Override
  public boolean integrity() {
    return Exam.integrity(mEnv);
  }
}
