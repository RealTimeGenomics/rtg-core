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
import com.rtg.util.Utils;

/**
 * Maintains a set of counts over the possible nucleotides.
 * <br/><br/>
 * The methods should be called in the order:<br/>
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
   * alleles we have counts for
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
    for (int i = 0; i < mDescription.size(); i++) {
      final int c = count(i);
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
    for (int i = 0; i < mDescription.size(); i++) {
      sb.append(" [").append(i).append("]  ");
      final int c = count(i);
      final String fc = Integer.toString(c);
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
  public abstract int count(final int index);

  /**
   * Get the current accumulated for the specified index.
   * @param index whose value to get.
   * @return the accumulated error.
   */
  public abstract double error(final int index);

  abstract Double alleleBalanceHomozygous(int allele, int total);

  abstract Double alleleBalance(int allele1, int allele2);

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
