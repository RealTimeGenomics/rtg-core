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

import java.util.Arrays;

import com.rtg.util.machine.MachineType;
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
   * @param clippedStart start position on read when excluding soft clipped bases (0 based)
   * @param clippedEnd end position on read when excluding soft clipped bases (0 based exclusive)
   */
  public AlignmentEnvironmentRead(VariantAlignmentRecord sam, VariantParams params, MachineType me, int clippedStart, int clippedEnd) {
    this(sam.getStart(), alignmentRecordToRead(sam, clippedStart, clippedEnd), alignmentRecordToQuality(sam, params, clippedStart, clippedEnd));
    assert me == null || !me.isCG();
    //SamUtils.checkCG(sam); // todo
  }

  /**
   * @param sam the SAM record
   * @param params the variant params object
   * @param me machine errors.
   */
  public AlignmentEnvironmentRead(VariantAlignmentRecord sam, VariantParams params, MachineType me) {
    this(sam.getStart(), alignmentRecordToRead(sam), alignmentRecordToQuality(sam, params));
  }

  private AlignmentEnvironmentRead(int start, byte[] read, double[] quality) {
    super(start);
    mRead = read;
    mQuality = quality;
  }

  private static byte[] alignmentRecordToRead(VariantAlignmentRecord sam) {
    return sam.getRead();
  }
  private static byte[] alignmentRecordToRead(VariantAlignmentRecord sam, int clippedStart, int clippedEnd) {
    return Arrays.copyOfRange(sam.getRead(), clippedStart, clippedEnd); // DnaUtils.encodeArrayCopy(sam.getRead(), clippedStart, clippedEnd - clippedStart); // DNA.byteDNAtoByte(sam.getRead());
  }

  private static double[] alignmentRecordToQuality(VariantAlignmentRecord sam, VariantParams params, int clippedStart, int clippedEnd) {
    final double[] ret = new double[clippedEnd - clippedStart];
    final byte[] quality = sam.getRecalibratedQuality();
    if (quality.length == 0) {
      final double qDef = VariantUtils.phredToProb(params.qDefault());
      for (int i = 0; i < ret.length; ++i) {
        ret[i] = qDef;
      }
    } else {
      for (int i = 0; i < ret.length; ++i) {
        ret[i] = VariantUtils.phredToProb(quality[i + clippedStart]);
      }
    }
    return ret;
  }

  private static double[] alignmentRecordToQuality(VariantAlignmentRecord sam, VariantParams params) {
    return alignmentRecordToQuality(sam, params, 0, sam.getRead().length);
  }

  @Override
  public int templateLength() {
    throw new UnsupportedOperationException();
  }

  @Override
  public boolean isInverted() {
    return false;
  }
}
