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

import junit.framework.TestCase;

/**
 */
public class GenotypeQualityOverDepthAnnotationTest extends TestCase {

  public void testName() {
    final GenotypeQualityOverDepthAnnotation gqdAnn = new GenotypeQualityOverDepthAnnotation();
    assertEquals("GQD", gqdAnn.getName());
    assertEquals("GQ / DP for a single sample", gqdAnn.getDescription());
    assertEquals(AnnotationDataType.DOUBLE, gqdAnn.getType());
  }

  public void test() {
    final GenotypeQualityOverDepthAnnotation gqdAnn = new GenotypeQualityOverDepthAnnotation();
    final VcfRecord rec = new VcfRecord();
    rec.setNumberOfSamples(6);
    rec.addFormatAndSample("GQ", "30");
    rec.addFormatAndSample("GQ", "40");
    rec.addFormatAndSample("GQ", "50");
    rec.addFormatAndSample("GQ", ".");
    rec.addFormatAndSample("GQ", "10");
    rec.addFormatAndSample("GQ", "10");
    rec.addFormatAndSample("DP", "50");
    rec.addFormatAndSample("DP", "40");
    rec.addFormatAndSample("DP", "30");
    rec.addFormatAndSample("DP", "10");
    rec.addFormatAndSample("DP", ".");
    rec.addFormatAndSample("DP", "0");
    assertEquals(0.6, (double) gqdAnn.getValue(rec, 0), 0.0001);
    assertEquals(1.0, (double) gqdAnn.getValue(rec, 1), 0.0001);
    assertEquals(1.6666, (double) gqdAnn.getValue(rec, 2), 0.0001);
    assertNull(gqdAnn.getValue(rec, 3));
    assertNull(gqdAnn.getValue(rec, 4));
    assertEquals(Double.POSITIVE_INFINITY, (double) gqdAnn.getValue(rec, 5), 0.01);
    assertNull(gqdAnn.getValue(rec, 6));
    rec.setNumberOfSamples(7);
    rec.addFormatAndSample("GQ", "11");
    assertNull(gqdAnn.getValue(rec, 0));
    assertEquals("Derived annotation GQD missing required fields in VCF header (FORMAT fields: DP GQ)", gqdAnn.checkHeader(null));
  }
}
