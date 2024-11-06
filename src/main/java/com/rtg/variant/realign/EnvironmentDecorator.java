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
