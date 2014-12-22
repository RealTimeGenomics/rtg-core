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

package com.rtg.vcf.annotation;

import java.util.List;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.VcfUtils;
import com.rtg.vcf.header.VcfHeader;

/**
 * Common code for inbreeding coefficient and Hardy Weinberg equilibrium probability
 * @see <a href="http://en.wikipedia.org/wiki/Hardy-Weinberg_principle#Inbreeding_coefficient">Inbreeding coefficient</a>
 */
@TestClass(value = {"com.rtg.vcf.annotation.InbreedingCoefficientAnnotationTest", "com.rtg.vcf.annotation.EquilibriumProbabilityAnnotationTest"})
public abstract class AbstractInbreedingCoefficientAnnotation extends AbstractDerivedAnnotation {

  /**
   * @param name the name for this annotation
   * @param description the description for this annotation
   */
  public AbstractInbreedingCoefficientAnnotation(String name, String description) {
    super(name, description, AnnotationDataType.DOUBLE);
  }

  @Override
  public Object getValue(VcfRecord record, int sampleNumber) {
    final List<String> gtList = record.getFormatAndSample().get(VcfUtils.FORMAT_GENOTYPE);
    if (gtList == null) {
      return null;
    }
    final int numAlleles = record.getAltCalls().size() + 1; //+1 for Reference allele
    final int[] alleleFreqCount = new int[numAlleles];
    int hetCount = 0;
    int total = 0;
    boolean diploid = false;
    for (String aGtList : gtList) {
      final int[] gts = VcfUtils.splitGt(aGtList);
      if (gts.length == 2) {
        diploid = true;
        if (gts[0] != gts[1]) {
          hetCount++;
        }
        for (int gt : gts) {
          if (gt < 0) {
            gt = 0;
          }
          alleleFreqCount[gt]++;
        }
        total++;
      }
    }
    if (!diploid) {
      return null; //Undefined for Haploid
    }
    return getValue(total, hetCount, getExpectedHetProb(total, alleleFreqCount));
  }

  protected double getExpectedHetProb(int total, int... haploidAlleleCount) {
    final int numAlleles = haploidAlleleCount.length;
    final double[] alleleFreqs = new double[numAlleles];
    for (int i = 0; i < numAlleles; i++) {
      alleleFreqs[i] = haploidAlleleCount[i] / (2.0 * total);
    }
    double expectedHetProb = 0;
    for (int i = 0; i < numAlleles; i++) {
      for (int j = i + 1; j < numAlleles; j++) {
        expectedHetProb += 2 * alleleFreqs[i] * alleleFreqs[j];
      }
    }
    return expectedHetProb;
  }


  @Override
  public String checkHeader(VcfHeader header) {
    return checkHeader(header, null, new String[]{VcfUtils.FORMAT_GENOTYPE});
  }

  protected abstract Double getValue(int total, int hetCount, double expectedHetProbability);

}
