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

package com.rtg.simulation.variants;

import com.rtg.mode.DNA;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Stores the specifics of a generated diploid mutation
 */
public class MutatorResult extends IntegralAbstract {
  private static final char[] DNA_CHARS = DNA.valueChars();
  private final byte[] mFirstHaplotype;
  private final byte[] mSecondHaplotype;
  private final int mConsumed;

  /**
   * @param firstHaplotype the sequence generated for the first haplotype
   * @param secondHaplotype the sequence generated for the second haplotype
   * @param consumed how many bases on the template were consumed
   */
  MutatorResult(byte[] firstHaplotype, byte[] secondHaplotype, int consumed) {
    super();
    mFirstHaplotype = firstHaplotype.clone();
    mSecondHaplotype = secondHaplotype.clone();
    mConsumed = consumed;
    assert integrity();
  }

  /**
   * @return Returns the first haplotype.
   */
  public byte[] getFirstHaplotype() {
    return mFirstHaplotype.clone();
  }

  /**
   * @return Returns the second haplotype.
   */
  public byte[] getSecondHaplotype() {
    return mSecondHaplotype.clone();
  }
  /**
   * Get consumed.
   * @return Returns the consumed.
   */
  public int getConsumed() {
    return mConsumed;
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mConsumed >= 0);
    Exam.assertNotNull(mFirstHaplotype);
    Exam.assertNotNull(mSecondHaplotype);
    checkHaplotype(mFirstHaplotype);
    checkHaplotype(mSecondHaplotype);
    return true;
  }

  static void checkHaplotype(byte[] nt) {
    for (final byte b : nt) {
      Exam.assertTrue(0 <= b && b <= 4);
    }
  }

  static void haplotypeToString(final StringBuilder sb, final byte[] haplotypes) {
    for (byte haplotype : haplotypes) {
      sb.append(DNA_CHARS[haplotype]);
    }
  }

  @Override
  public void toString(StringBuilder sb) {
    sb.append(mConsumed).append(":");
    haplotypeToString(sb, mFirstHaplotype);
    sb.append(":");
    haplotypeToString(sb, mSecondHaplotype);
  }

}
