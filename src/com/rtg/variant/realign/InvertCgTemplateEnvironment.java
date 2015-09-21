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
