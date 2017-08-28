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
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.machine.MachineType;
import com.rtg.variant.VariantAlignmentRecord;
import com.rtg.variant.VariantParams;
import com.rtg.variant.util.CgUnroller;
import com.rtg.variant.util.VariantUtils;

/**
 * Wrap a CG alignment record so that prepared for use in an environment
 * for all-paths calculations.
 */
public class AlignmentEnvironmentCG extends AbstractAlignmentEnvironment {

  private final boolean mIsInverted;

  /**
   * @param var the alignment record
   * @param params the variant params object
   * @param template reference nucleotides
   * @param machineType which machine type
   */
  public AlignmentEnvironmentCG(final VariantAlignmentRecord var, final VariantParams params, final byte[] template, MachineType machineType) {
    super(var.getStart());
    assert machineType.isCG();
    final CgUnroller.OrientedRead orient = CgUnroller.unrollCgRead(var, template);
    if (orient == null) {
      throw new NoTalkbackSlimException("Invalid CG alignment. Could not reconstruct original read. record=" + var);
    }
    mRead = orient.getRead();
    final int len = mRead.length;
    //System.err.println("CG len=" + len);
    if (machineType == MachineType.COMPLETE_GENOMICS && len != CgUtils.CG_RAW_READ_LENGTH) {
      throw new NoTalkbackSlimException("Invalid CG version 1 alignment. Unexpected reconstructed read length. record=" + var);
    }
    if (machineType == MachineType.COMPLETE_GENOMICS_2 && len != CgUtils.CG2_RAW_READ_LENGTH) {
      throw new NoTalkbackSlimException("Invalid CG version 2 alignment. Unexpected reconstructed read length. record=" + var);
    }
    mQuality = new double[len];
    final byte[] qChar = orient.getQuality();
    if (qChar == null) {
      final double qDef = VariantUtils.phredToProb(params.qDefault());
      for (int i = 0; i < len; ++i) {
        mQuality[i] = qDef;
      }
    } else {
      for (int i = 0; i < len; ++i) {
        mQuality[i] = VariantUtils.phredToProb(qChar[i]);
      }
    }
    mIsInverted = orient.isInverted();
  }

  @Override
  public int templateLength() {
    throw new UnsupportedOperationException();
  }

  @Override
  public boolean isInverted() {
    return mIsInverted;
  }
}
