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

import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.VcfUtils;

import junit.framework.TestCase;

/**
 */
public class EquilibriumProbabilityAnnotationTest extends TestCase {

  public void testName() {
    final EquilibriumProbabilityAnnotation epAnn = new EquilibriumProbabilityAnnotation();
    assertEquals("EP", epAnn.getName());
    assertEquals("Phred scaled probability that site is not in Hardy-Weinberg equilibrium", epAnn.getDescription());
    assertEquals(AnnotationDataType.DOUBLE, epAnn.getType());
    assertEquals("Derived annotation EP missing required fields in VCF header (FORMAT fields: GT)", epAnn.checkHeader(null));
  }

  public void testHaploidNoCalculation() {
    final EquilibriumProbabilityAnnotation epAnn = new EquilibriumProbabilityAnnotation();
    final VcfRecord rec = new VcfRecord();
    rec.setRefCall("A")
        .addAltCall("C")
        .setNumberOfSamples(3)
        .addFormatAndSample(VcfUtils.FORMAT_GENOTYPE, "0")
        .addFormatAndSample(VcfUtils.FORMAT_GENOTYPE, ".")
        .addFormatAndSample(VcfUtils.FORMAT_GENOTYPE, "1");
    assertNull(epAnn.getValue(rec, -2));
  }

  public void testCoefficientCalculation() {
    final EquilibriumProbabilityAnnotation epAnn = new EquilibriumProbabilityAnnotation();
    VcfRecord rec = InbreedingCoefficientAnnotationTest.makeTwoAlleleRecord(1469, 138, 5);
    assertEquals(0.0554, (Double) epAnn.getValue(rec, 0), 0.001);
    rec = InbreedingCoefficientAnnotationTest.makeTwoAlleleRecord(0, 11, 0);
    assertEquals(23.886, (Double) epAnn.getValue(rec, 0), 0.001);
    rec = InbreedingCoefficientAnnotationTest.makeTwoAlleleRecord(1, 0, 1);
    assertEquals(4.343, (Double) epAnn.getValue(rec, 0), 0.001);
    rec = InbreedingCoefficientAnnotationTest.makeTwoAlleleRecord(5, 11, 5);
    assertEquals(0.103, (Double) epAnn.getValue(rec, 0), 0.001);
    rec = InbreedingCoefficientAnnotationTest.makeTwoAlleleRecord(1, 1, 1);
    assertEquals(0.724, (Double) epAnn.getValue(rec, 0), 0.001);
    rec = InbreedingCoefficientAnnotationTest.makeTwoAlleleRecord(0, 0, 10);
    assertEquals(0.0, (Double) epAnn.getValue(rec, 0), 0.001);
    rec = InbreedingCoefficientAnnotationTest.makeTwoAlleleRecord(10, 0, 0);
    assertEquals(0.0, (Double) epAnn.getValue(rec, 0), 0.001);
    rec = InbreedingCoefficientAnnotationTest.makeThreeAlleleRecord(1208, 222, 1146, 10, 110, 288);
    assertEquals(0.639, (Double) epAnn.getValue(rec, -1), 0.001);
    rec = InbreedingCoefficientAnnotationTest.makeThreeAlleleRecord(0, 0, 0, 1469, 138, 5);
    assertEquals(0.0554, (Double) epAnn.getValue(rec, 123), 0.001);
  }
}
