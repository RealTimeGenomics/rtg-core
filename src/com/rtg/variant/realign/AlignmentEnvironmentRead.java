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
    //XXX maybe check whether we can pass this directly instead of copying.
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
