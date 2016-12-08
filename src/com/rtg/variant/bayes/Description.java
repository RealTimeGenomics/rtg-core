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

package com.rtg.variant.bayes;

import com.reeltwo.jumble.annotations.TestClass;

/**
 * Provides names for a set of haploid hypotheses. e.g. alleles or haplotypes
 */
@TestClass("com.rtg.variant.bayes.snp.DescriptionSnpTest")
public abstract class Description {

  /**
   * Name of the hypothesis selected by index.
   * @param index select hypothesis (0 based).
   * @return the name.
   */
  public abstract String name(int index);

  /**
   * Return the index of the given name
   * @param name the name to search for
   * @return the index of name, or -1 if the name cannot be found
   */
  public int indexOf(String name) {
    for (int i = 0; i < size(); ++i) {
      if (name(i).equals(name)) {
        return i;
      }
    }
    return -1;
  }

  /**
   * Write the name of the hypothesis as a sequence of bytes.
   * @param sb where the name is to be written.
   * @param index select hypothesis (0 based).
   */
  public void writeName(StringBuilder sb, int index) {
    if (name(index) != null) {
      sb.append(name(index));
    }
  }

  /**
   * @return the number of different hypotheses.
   */
  public abstract int size();

  /**
   * Check if the index is in valid range.
   * @param index to be checked.
   * @return true iff index is in range.
   */
  public boolean valid(int index) {
    return 0 <= index && index < size();
  }

  /**
   * @return the minimum length (in nt) of any of the haploid hypotheses. <b>Integer.MAX_VALUE</b> if hypotheses is empty.
   */
  public abstract int minLength();

  /**
   * @return the maximum length (in nt) of any of the haploid hypotheses. 0 if hypotheses is empty.
   */
  public abstract int maxLength();

}
