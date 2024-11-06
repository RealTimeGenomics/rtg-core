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
