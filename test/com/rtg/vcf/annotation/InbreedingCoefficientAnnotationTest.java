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
 * Produces an info field containing the inbreeding coefficient for a set of samples.
 * @see <a href="http://en.wikipedia.org/wiki/Hardy-Weinberg_principle#Inbreeding_coefficient">Inbreeding coefficient</a>
 */
public class InbreedingCoefficientAnnotationTest extends TestCase {

  public void testName() {
    final InbreedingCoefficientAnnotation icAnn = new InbreedingCoefficientAnnotation();
    assertEquals("IC", icAnn.getName());
    assertEquals("Inbreeding Coefficient", icAnn.getDescription());
    assertEquals(AnnotationDataType.DOUBLE, icAnn.getType());
    assertEquals("Derived annotation IC missing required fields in VCF header (FORMAT fields: GT)", icAnn.checkHeader(null));
  }

  public void testHaploidNoCalculation() {
    final InbreedingCoefficientAnnotation icAnn = new InbreedingCoefficientAnnotation();
    final VcfRecord rec = new VcfRecord();
    rec.setRefCall("A");
    rec.addAltCall("C");
    rec.setNumberOfSamples(3);
    rec.addFormatAndSample(VcfUtils.FORMAT_GENOTYPE, "0");
    rec.addFormatAndSample(VcfUtils.FORMAT_GENOTYPE, ".");
    rec.addFormatAndSample(VcfUtils.FORMAT_GENOTYPE, "1");
    assertNull(icAnn.getValue(rec, 0));
  }

  public void testCoefficientCalculation() {
    final InbreedingCoefficientAnnotation icAnn = new InbreedingCoefficientAnnotation();
    VcfRecord rec = makeTwoAlleleRecord(1469, 138, 5);
    assertEquals(0.023, (Double) icAnn.getValue(rec, 1), 0.001);
    rec = makeTwoAlleleRecord(0, 11, 0);
    assertEquals(-1.000, (Double) icAnn.getValue(rec, 23), 0.001);
    rec = makeTwoAlleleRecord(1, 0, 1);
    assertEquals(1.000, (Double) icAnn.getValue(rec, 0), 0.001);
    rec = makeTwoAlleleRecord(5, 11, 5);
    assertEquals(-0.048, (Double) icAnn.getValue(rec, -1), 0.001);
    rec = makeTwoAlleleRecord(1, 1, 1);
    assertEquals(0.333, (Double) icAnn.getValue(rec, 0), 0.001);
    rec = makeTwoAlleleRecord(0, 0, 10);
    assertNull(icAnn.getValue(rec, 0));
    rec = makeTwoAlleleRecord(10, 0, 0);
    assertNull(icAnn.getValue(rec, 0));
    rec = makeThreeAlleleRecord(1208, 222, 1146, 10, 110, 288);
    assertEquals(0.010, (Double) icAnn.getValue(rec, 0), 0.001);
    rec = makeThreeAlleleRecord(0, 0, 0, 1469, 138, 5);
    assertEquals(0.023, (Double) icAnn.getValue(rec, 0), 0.001);
    rec = makeOneAlleleRecord(11);
    assertNull(icAnn.getValue(rec, 22));
  }

  protected static VcfRecord makeOneAlleleRecord(int aa) {
    final VcfRecord rec = new VcfRecord();
    rec.setRefCall("A");
    rec.setNumberOfSamples(aa);
    for (int i = 0; i < aa / 2; i++) {
      rec.addFormatAndSample(VcfUtils.FORMAT_GENOTYPE, "0/0");
    }
    for (int i = aa / 2; i < aa; i++) {
      rec.addFormatAndSample(VcfUtils.FORMAT_GENOTYPE, "./.");
    }
    return rec;
  }

  protected static VcfRecord makeTwoAlleleRecord(int aa, int ac, int cc) {
    final VcfRecord rec = new VcfRecord();
    rec.setRefCall("A");
    rec.addAltCall("C");
    rec.setNumberOfSamples(aa + ac + cc);
    for (int i = 0; i < aa / 2; i++) {
      rec.addFormatAndSample(VcfUtils.FORMAT_GENOTYPE, "0/0");
    }
    for (int i = aa / 2; i < aa; i++) {
      rec.addFormatAndSample(VcfUtils.FORMAT_GENOTYPE, "./.");
    }
    for (int i = 0; i < ac / 2; i++) {
      rec.addFormatAndSample(VcfUtils.FORMAT_GENOTYPE, "0/1");
    }
    for (int i = ac / 2; i < ac; i++) {
      rec.addFormatAndSample(VcfUtils.FORMAT_GENOTYPE, "1/0");
    }
    for (int i = 0; i < cc; i++) {
      rec.addFormatAndSample(VcfUtils.FORMAT_GENOTYPE, "1/1");
    }
    return rec;
  }

  protected static VcfRecord makeThreeAlleleRecord(int aa, int ac, int ag, int cc, int cg, int gg) {
    final VcfRecord rec = new VcfRecord();
    rec.setRefCall("A");
    rec.addAltCall("C");
    rec.addAltCall("G");
    rec.setNumberOfSamples(aa + ac + ag + cc + cg + gg);
    for (int i = 0; i < aa; i++) {
      rec.addFormatAndSample(VcfUtils.FORMAT_GENOTYPE, "0/0");
    }
    for (int i = 0; i < ac / 2; i++) {
      rec.addFormatAndSample(VcfUtils.FORMAT_GENOTYPE, "0/1");
    }
    for (int i = ac / 2; i < ac; i++) {
      rec.addFormatAndSample(VcfUtils.FORMAT_GENOTYPE, "1/0");
    }
    for (int i = 0; i < ag / 2; i++) {
      rec.addFormatAndSample(VcfUtils.FORMAT_GENOTYPE, "0/2");
    }
    for (int i = ag / 2; i < ag; i++) {
      rec.addFormatAndSample(VcfUtils.FORMAT_GENOTYPE, "2/0");
    }
    for (int i = 0; i < cc; i++) {
      rec.addFormatAndSample(VcfUtils.FORMAT_GENOTYPE, "1/1");
    }
    for (int i = 0; i < cg / 2; i++) {
      rec.addFormatAndSample(VcfUtils.FORMAT_GENOTYPE, "1/2");
    }
    for (int i = cg / 2; i < cg; i++) {
      rec.addFormatAndSample(VcfUtils.FORMAT_GENOTYPE, "2/1");
    }
    for (int i = 0; i < gg; i++) {
      rec.addFormatAndSample(VcfUtils.FORMAT_GENOTYPE, "2/2");
    }
    return rec;
  }
}
