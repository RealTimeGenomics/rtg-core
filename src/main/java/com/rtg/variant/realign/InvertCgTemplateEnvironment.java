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

import com.rtg.reader.CgUtils;
import com.rtg.util.machine.MachineType;

/**
 * An environment with the read and the template reversed, used for CG
 * reads so the gap/overlap layout is consistent.
 *
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

  private final int mOffset;

  /**
   * @param env environment to be inverted.
   * @param mt the machine type, to determine whether to do inversion for version 1 or version 2 reads
   */
  public InvertCgTemplateEnvironment(final EnvironmentCombined env, MachineType mt) {
    super(env.mSamEnv, env.mSamEnv.start(), env.maxShift(), env.mTemEnv);
    if (mt == MachineType.COMPLETE_GENOMICS) {
      mOffset = CgUtils.CG_EXPECTED_LENGTH - 1;
    } else if (mt == MachineType.COMPLETE_GENOMICS_2) {
      mOffset = CgUtils.CG2_EXPECTED_LENGTH - 1;
    } else {
      throw new IllegalArgumentException("Unsupported machine type: " + mt);
    }
  }

  @Override
  public byte template(int index) {
    return super.template(mOffset - index);
  }

  @Override
  public int absoluteTemplatePosition(int index) {
    return super.absoluteTemplatePosition(mOffset - index);
  }


}
