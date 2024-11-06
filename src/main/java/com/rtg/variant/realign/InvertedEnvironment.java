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
