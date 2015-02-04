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

import com.rtg.mode.DNA;
import com.rtg.variant.AbstractMachineErrorParams;
import com.rtg.variant.VariantAlignmentRecord;
import com.rtg.variant.VariantParams;
import com.rtg.variant.util.VariantUtils;

/**
 * Wrap a SAM record so that prepared for use in an environment
 * for all-paths calculations.
 */
public final class AlignmentEnvironmentRead extends AbstractAlignmentEnvironment {

  /**
   * @param sam the SAM record
   * @param params the variant params object
   * @param me machine errors.
   */
  public AlignmentEnvironmentRead(VariantAlignmentRecord sam, VariantParams params, AbstractMachineErrorParams me) {
    super(sam.getStart());
    assert !me.isCG();
    //SamUtils.checkCG(sam); // todo
    mRead = DNA.byteDNAtoByte(sam.getRead());
    final int len = mRead.length;
    mQuality = new double[len];
    final byte[] quality = sam.getQuality();
    if (quality.length == 0) {
      final double qDef = VariantUtils.phredToProb(params.qDefault());
      for (int i = 0; i < len; i++) {
        mQuality[i] = qDef;
      }
    } else {
      for (int i = 0; i < len; i++) {
        final int phred = me.getPhred(quality[i], i);
        mQuality[i] = VariantUtils.phredToProb(phred);
      }
    }
  }

  @Override
  public int templateLength() {
    throw new UnsupportedOperationException();
  }

  @Override
  public boolean cgOverlapOnLeft() {
    return true;
  }
}
