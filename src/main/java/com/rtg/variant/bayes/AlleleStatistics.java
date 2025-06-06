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
import com.rtg.util.MathUtils;
import com.rtg.util.Utils;

/**
 * Maintains a set of counts over the possible nucleotides.
 * <br><br>
 * The methods should be called in the order:<br>
 * <code>((reset | constructor) increment* count* )*</code>
 */
@TestClass("com.rtg.variant.bayes.AlleleStatisticsIntTest")
public abstract class AlleleStatistics<T extends AlleleStatistics<T>> {

  private final Description mDescription;

  /**
   * Create a new empty count object.
   * @param description set of possible values
   */
  public AlleleStatistics(final Description description) {
    mDescription = description;
  }


  /**
   * Alleles we have counts for.
   * @return the description for the alleles
   */
  public final Description getDescription() {
    return mDescription;
  }

  /**
   * Write statistics output for use in variant output file.
   * @param sb string builder to write statistics to
   * @param separator between entries
   */
  public final void output(StringBuilder sb, char separator) {
    for (int i = 0; i < mDescription.size(); ++i) {
      final long c = MathUtils.round(count(i));
      if (c == 0) {
        continue;
      }
      sb.append(separator);
      final String name = mDescription.name(i);
      if (name != null) {
        sb.append(name);
      }
      sb.append(separator);
      sb.append(c);
      sb.append(separator);
      final double qe = error(i);
      sb.append(Utils.realFormat(qe, 3));
    }
  }

  @Override
  public String toString() {
    final StringBuilder sb = new StringBuilder();
    for (int i = 0; i < mDescription.size(); ++i) {
      sb.append(" [").append(i).append("]  ");
      final long c = MathUtils.round(count(i));
      final String fc = Long.toString(c);
      sb.append(fc);

      sb.append("  ");
      final double qe = error(i);
      final String fqe = String.format("%1$04.3f", qe);
      sb.append(fqe);
    }
    return sb.toString();
  }

  /**
   * Get the current count for the specified index.
   * @param index whose value to get.
   * @return the count.
   */
  public abstract double count(final int index);

  /**
   * Get the current forward count for the specified index.
   * @param index whose value to get.
   * @return the count.
   */
  public double forward(final int index) {
    return forward1(index) + forward2(index);
  }

  /**
   * Get the current R1 forward count for the specified index.
   * @param index whose value to get.
   * @return the count.
   */
  public abstract double forward1(final int index);

  /**
   * Get the current R2 forward count for the specified index.
   * @param index whose value to get.
   * @return the count.
   */
  public abstract double forward2(final int index);

  /**
   * Get the current backward count for the specified index.
   * @param index whose value to get.
   * @return the count.
   */
  public double backward(final int index) {
    return backward1(index) + backward2(index);
  }

  /**
   * Get the current R1 backward count for the specified index.
   * @param index whose value to get.
   * @return the count.
   */
  public abstract double backward1(final int index);

  /**
   * Get the current R2 backward count for the specified index.
   * @param index whose value to get.
   * @return the count.
   */
  public abstract double backward2(final int index);

  /**
   * Get the current accumulated for the specified index.
   * @param index whose value to get.
   * @return the accumulated error.
   */
  public abstract double error(final int index);

  /**
   * Get the current accumulated for the specified index.
   * @param index whose value to get.
   * @return the accumulated error.
   */
  public abstract double qa(final int index);


  Double alleleBalanceHomozygous(int allele, int total) {
    return MathUtils.hoeffdingPhred(total, count(allele), 1.0);
  }

  Double alleleBalance(int allele1, int allele2) {
    final double trials = count(allele1) + count(allele2);
    final double observed = count(allele1);
    return MathUtils.hoeffdingPhred(trials, observed, 0.5);
  }

  abstract Double strandBias(int allele);

  abstract Double unmatedBias(int allele, double unmatedProbability);

  /**
   * Returns a new AlleleStatistics that has copied counts from the current
   * description into a new description according to specified description mappings.
   * @param newDescription the new description
   * @param mapping a mapping array. counts for old description <code>i</code> should
   * map to new description <code>mapping[i]</code>
   * @return T new AlleleStatistics instance with remapped entries
   */
  public abstract T remap(Description newDescription, int[] mapping);
}
