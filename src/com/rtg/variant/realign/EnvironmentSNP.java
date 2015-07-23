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

/**
 * Adds a SNP at one position along the template.
 * Limitation: The <code>toString</code> method is not updated,
 * so still shows the underlying template.
 *
 */
public class EnvironmentSNP extends EnvironmentDecorator {

  private final int mPosition;
  private final byte mNucleotide;

  /**
   * @param env the original environment
   * @param tpos the relative template position to change (0 = start of read)
   * @param nucleotide the new template nucleotide at that position.
   */
  EnvironmentSNP(Environment env, int tpos, byte nucleotide) {
    super(env);
    mPosition = tpos;
    mNucleotide = nucleotide;
  }

  @Override
  public byte template(int index) {
    if (index == mPosition) {
      return mNucleotide;
    } else {
      return super.template(index);
    }
  }
}
