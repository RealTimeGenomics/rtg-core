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
public class PloidyAnnotationTest extends TestCase {

  public void testName() {
    final PloidyAnnotation pdAnn = new PloidyAnnotation();
    assertEquals("PD", pdAnn.getName());
    assertEquals("Ploidy of sample. 'h'=>haploid, 'd'=>diploid", pdAnn.getDescription());
    assertEquals(AnnotationDataType.STRING, pdAnn.getType());
  }

  public void test() {
    final PloidyAnnotation pdAnn = new PloidyAnnotation();
    final VcfRecord rec = new VcfRecord();
    rec.setNumberOfSamples(6);
    rec.addFormatAndSample("GT", "0/1");
    rec.addFormatAndSample("GT", ".");
    rec.addFormatAndSample("GT", "1/1");
    rec.addFormatAndSample("GT", "./.");
    rec.addFormatAndSample("GT", "1");
    rec.addFormatAndSample("GT", "0");
    assertEquals("d", pdAnn.getValue(rec, 0));
    assertNull(pdAnn.getValue(rec, 1));
    assertEquals("d", pdAnn.getValue(rec, 2));
    assertNull(pdAnn.getValue(rec, 3));
    assertEquals("h", pdAnn.getValue(rec, 4));
    assertEquals("h", pdAnn.getValue(rec, 5));
    assertNull(pdAnn.getValue(rec, 6));
    assertEquals("Derived annotation PD missing required fields in VCF header (FORMAT fields: GT)", pdAnn.checkHeader(null));
  }
}
