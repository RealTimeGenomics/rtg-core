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
import com.rtg.sam.SamUtils;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.variant.AbstractMachineErrorParams;
import com.rtg.variant.VariantAlignmentRecord;
import com.rtg.variant.VariantParams;
import com.rtg.variant.util.CgUnroller;
import com.rtg.variant.util.VariantUtils;

/**
 * Wrap a CG alignment record so that prepared for use in an environment
 * for all-paths calculations.
 */
public class AlignmentEnvironmentCG extends AbstractAlignmentEnvironment {

  private final boolean mCgOverlapOnLeft;

  /**
   * @param var the alignment record
   * @param params the variant params object
   * @param me machine errors.
   * @param template reference nucleotides
   */
  public AlignmentEnvironmentCG(final VariantAlignmentRecord var, final VariantParams params, final AbstractMachineErrorParams me, final byte[] template) {
    super(var.getStart());
    assert me.isCG();
    final CgUnroller.OrientedRead orient = CgUnroller.unrollCgRead(var, template);
    if (orient == null) {
      throw new NoTalkbackSlimException("Invalid CG alignment record=" + var.toString());
    }
    mRead = DNA.byteDNAtoByte(orient.getRead());
    final int len = mRead.length;
    //System.err.println("CG len=" + len);
    if (len != SamUtils.CG_RAW_READ_LENGTH) {
      throw new NoTalkbackSlimException("Invalid CG alignment record=" + var.toString());
    }
    mQuality = new double[len];
    final byte[] qChar = orient.getQuality();
    final String quality = qChar == null ? null : new String(qChar);
    if (quality == null) {
      final double qDef = VariantUtils.phredToProb(params.qDefault());
      for (int i = 0; i < len; i++) {
        mQuality[i] = qDef;
      }
    } else {
      for (int i = 0; i < len; i++) {
        final char ch = quality.charAt(i);
        final int phred = me.getPhred(ch, i);
        mQuality[i] = VariantUtils.phredToProb(phred);
      }
    }
    mCgOverlapOnLeft = orient.isCgOverlapOnLeft();
  }

  @Override
  public int templateLength() {
    throw new UnsupportedOperationException();
  }

  @Override
  public boolean cgOverlapOnLeft() {
    return mCgOverlapOnLeft;
  }
}
