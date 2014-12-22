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

import com.rtg.util.TestUtils;

import junit.framework.TestCase;

/**
 */
public class DerivedAnnotationsTest extends TestCase {

  public void testEnum() {
    TestUtils.testEnum(DerivedAnnotations.class, "[IC, EP, LAL, QD, NAA, AC, AN, GQD, ZY, PD]");
    assertTrue(DerivedAnnotations.IC.getAnnotation() instanceof InbreedingCoefficientAnnotation);
    assertTrue(DerivedAnnotations.EP.getAnnotation() instanceof EquilibriumProbabilityAnnotation);
    assertTrue(DerivedAnnotations.LAL.getAnnotation() instanceof LongestAlleleAnnotation);
    assertTrue(DerivedAnnotations.QD.getAnnotation() instanceof QualOverDepthAnnotation);
    assertTrue(DerivedAnnotations.NAA.getAnnotation() instanceof NumberOfAltAllelesAnnotation);
    assertTrue(DerivedAnnotations.AC.getAnnotation() instanceof AlleleCountInGenotypesAnnotation);
    assertTrue(DerivedAnnotations.AN.getAnnotation() instanceof NumberAllelesInGenotypesAnnotation);
    assertTrue(DerivedAnnotations.GQD.getAnnotation() instanceof GenotypeQualityOverDepthAnnotation);
    assertTrue(DerivedAnnotations.ZY.getAnnotation() instanceof ZygosityAnnotation);
    assertTrue(DerivedAnnotations.PD.getAnnotation() instanceof PloidyAnnotation);
  }

}
