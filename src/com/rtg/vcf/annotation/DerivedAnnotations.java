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

import java.util.EnumSet;

/**
 * Enumeration of the derived annotations.
 */
public enum DerivedAnnotations {
  /** Inbreeding Coefficient */
  IC(new InbreedingCoefficientAnnotation()),
  /** Hardy-Weinberg Equilibrium Probability */
  EP(new EquilibriumProbabilityAnnotation()),
  /** Length of the longest allele */
  LAL(new LongestAlleleAnnotation()),
  /** QUAL / DP */
  QD(new QualOverDepthAnnotation()),
  /** Number of alternative alleles */
  NAA(new NumberOfAltAllelesAnnotation()),
  /**
   * Allele count in genotypes, for each alternative allele
   * Note: Is a multiple value annotation (can not be used for AVR ML)
   */
  AC(new AlleleCountInGenotypesAnnotation()),
  /** Total number of alleles in called genotypes */
  AN(new NumberAllelesInGenotypesAnnotation()),
  /** GQ / DP for a single sample */
  GQD(new GenotypeQualityOverDepthAnnotation()),
  /** Zygosity */
  ZY(new ZygosityAnnotation()),
  /** Ploidy */
  PD(new PloidyAnnotation());

  private final AbstractDerivedAnnotation mAnnotation;

  private DerivedAnnotations(AbstractDerivedAnnotation annotation) {
    assert name().equals(annotation.getName());
    mAnnotation = annotation;
  }

  /**
   * Get the annotation associated with this enum value.
   * @return the annotation associated with this enum value.
   */
  public AbstractDerivedAnnotation getAnnotation() {
    return mAnnotation;
  }

  /**
   * Get the set of derived annotations that produce a single value.
   * @return the set of derived annotations that produce a single value.
   */
  public static EnumSet<DerivedAnnotations> singleValueAnnotations() {
    return EnumSet.complementOf(EnumSet.of(AC));
  }

  /**
   * Get the set of derived annotations that produce a single numeric value.
   * @return the set of derived annotations that produce a single numeric value.
   */
  public static EnumSet<DerivedAnnotations> singleValueNumericAnnotations() {
    final EnumSet<DerivedAnnotations> ret = EnumSet.noneOf(DerivedAnnotations.class);
    for (final DerivedAnnotations ann : singleValueAnnotations()) {
      if (ann.getAnnotation().getType() == AnnotationDataType.INTEGER || ann.getAnnotation().getType() == AnnotationDataType.DOUBLE) {
        ret.add(ann);
      }
    }
    return ret;
  }
}
