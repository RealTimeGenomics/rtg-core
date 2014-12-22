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
public class ZygosityAnnotationTest extends TestCase {

  public void testName() {
    final ZygosityAnnotation zyAnn = new ZygosityAnnotation();
    assertEquals("ZY", zyAnn.getName());
    assertEquals("Zygosity of sample. 'e'=>heterozygous, 'o'=>homozygous", zyAnn.getDescription());
    assertEquals(AnnotationDataType.STRING, zyAnn.getType());
  }

  public void test() {
    final ZygosityAnnotation zyAnn = new ZygosityAnnotation();
    final VcfRecord rec = new VcfRecord();
    rec.setNumberOfSamples(3);
    rec.addFormatAndSample("GT", "0/1");
    rec.addFormatAndSample("GT", ".");
    rec.addFormatAndSample("GT", "1/1");
    assertEquals("e", zyAnn.getValue(rec, 0));
    assertNull(zyAnn.getValue(rec, 1));
    assertEquals("o", zyAnn.getValue(rec, 2));
    assertNull(zyAnn.getValue(rec, 3));
    assertEquals("Derived annotation ZY missing required fields in VCF header (FORMAT fields: GT)", zyAnn.checkHeader(null));
  }
}
