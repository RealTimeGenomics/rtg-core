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
package com.rtg.variant;

import java.io.IOException;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.reader.FastaUtils;
import com.rtg.util.Params;
import com.rtg.util.machine.MachineType;
import com.rtg.variant.realign.RealignParams;

/**
 * Abstract machine error parameters class
 */
@TestClass("com.rtg.variant.MachineErrorParamsTest")
public abstract class AbstractMachineErrorParams implements Params {

  /**
   * Get a phred score from an ASCII quality character optionally
   * correcting it.
   * @param qualChar  original quality character.
   * @param readPos position on read of <code>qualChar</code>
   * @return the possibly corrected phred score.
   */
  public final int getPhred(final char qualChar, int readPos) {
    assert qualChar >= FastaUtils.PHRED_LOWER_LIMIT_CHAR;
    final byte rawQuality = (byte) (qualChar - FastaUtils.PHRED_LOWER_LIMIT_CHAR);
    return getPhred(rawQuality, readPos);
  }

  /**
   * Get a phred score from a binary quality value optionally
   * correcting it.
   * @param quality original quality value.
   * @param readPos position on read of <code>qualChar</code>
   * @return the possibly corrected phred score.
   */
  public abstract int getPhred(final byte quality, int readPos);

  /**
   * Get the flag to indicate if CG outer base trimming is to be used.
   * @return the flag - true if trimming is to be used.
   */
  public boolean cgTrimOuterBases() {
    return machineType() == MachineType.COMPLETE_GENOMICS && MachineErrorParamsBuilder.CG_TRIM;
  }

  /**
   * Get the CG v1 small gap distribution for 0,1,2,3.
   * @return the gap distribution, always non-null.
   */
  public abstract double[] smallGapDistribution();

  /**
   * Get the CG v1 large gap distribution for 4,5,6,7,8.
   * @return the gap distribution, always non-null.
   */
  public abstract double[] gapDistribution();

  /**
   * Get the CG v1 overlap probability distribution for -4,-3,-2,-1,0.
   *
   * @return the overlap distribution, always non-null.
   */
  public abstract double[] overlapDistribution();

  /**
   * Get the CG v2 overlap probability distribution for -7,-6,-5,-4,-3,-2,-1,0.
   *
   * @return the overlap distribution, always non-null.
   */
  public abstract double[] overlapDistribution2();

  /**
   * Return the machine type specified in the priors
   *
   * @return the expected machine type
   */
  public abstract MachineType machineType();

  /**
   * True if the CG priors (in particular the overlap distribution) has been explicitly set.
   *
   * @return true for CG.
   */
  public abstract boolean isCG();

  /**
   * Get the optional quality calibration curve.
   *
   * @return null, or an array with 64 entries, each in the range 0..63.
   */
  protected abstract int[] qualityCurve();

  /**
   * Get the length distribution of sequencing machine deletion errors.
   * This does not include the underlying mutation distribution.
   * So the total deletion distribution should be the weighted sum
   * of this one and <code>deletionDistribution()</code>.
   *
   * @return an array that sums to 1.0.  Entry 0 will be 0.0.
   */
  public abstract double[] errorDelDistribution();

  /**
   * Get the sequencing machine error event rate for deletion events.
   * @return a probability between 0 and 1.
   */
  public abstract double errorDelEventRate();

  /**
   * Get the sequencing machine error rate for bases deleted.
   * @return a probability between 0 and 1.
   */
  public abstract double errorDelBaseRate();

  /**
   * Get the length distribution of sequencing machine insertion errors.
   * This does not include the underlying mutation distribution.
   * So the total insertion distribution should be the weighted sum
   * of this one and <code>insertionDistribution()</code>.
   *
   * @return an array that sums to 1.0.  Entry 0 will be 0.0.
   */
  public abstract double[] errorInsDistribution();

  /**
   * Get the sequencing machine error rate for insert events.
   * @return a probability between 0 and 1.
   */
  public abstract double errorInsEventRate();

  /**
   * Get the sequencing machine error rate for bases inserted.
   * @return a probability between 0 and 1.
   */
  public abstract double errorInsBaseRate();

  /**
   * Get the length distribution of sequencing machine MNP errors.
   * This includes SNP errors (when length equals 1).
   * This does not include the underlying mutation distribution.
   * So the total MNP distribution should be the weighted sum
   * of this one and <code>insertionDistribution()</code>.
   *
   * @return an array that sums to 1.0.  Entry 0 will be 0.0.
   */
  public abstract double[] errorMnpDistribution();

  /**
   * Get the sequencing machine error rate for an MNP event starting.
   * Since an MNP of length 1 is a SNP, this includes the probability
   * of a SNP error.  In fact, the SNP rate is roughly equal to this
   * MNP event rate times the average length of an MNP.
   *
   * @return a probability between 0 and 1.
   */
  public abstract double errorMnpEventRate();

  /**
   * Get the sequencing machine error rate for single bases changed.
   * This is calculated as the MNP event rate * the average length of a MNP.
   * @return a probability between 0 and 1.
   */
  public abstract double errorSnpRate();

  /**
   * @return the RealignParams associated with these errors.
   */
  public abstract RealignParams realignParams();

  @Override
  public void close() throws IOException { }

}
