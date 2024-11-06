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
