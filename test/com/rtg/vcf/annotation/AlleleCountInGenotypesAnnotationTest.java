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
public class AlleleCountInGenotypesAnnotationTest extends TestCase {

  public void testName() {
    final AlleleCountInGenotypesAnnotation ann = new AlleleCountInGenotypesAnnotation();
    assertEquals("AC", ann.getName());
    assertEquals("Allele count in genotypes, for each alternative allele, in the same order as listed", ann.getDescription());
    assertEquals(AnnotationDataType.INTEGER, ann.getType());
  }

  public void test() {
    final AlleleCountInGenotypesAnnotation ann = new AlleleCountInGenotypesAnnotation();
    VcfRecord record = new VcfRecord();
    record.setNumberOfSamples(4);
    record.setRefCall("G");
    record.addAltCall("A");
    record.addAltCall("C");
    record.addFormatAndSample("GT", "1/0");
    record.addFormatAndSample("GT", "1/1");
    record.addFormatAndSample("GT", ".");
    record.addFormatAndSample("GT", "2");
    final int[] vals = (int[]) ann.getValue(record, -1);
    assertNotNull(vals);
    assertEquals(2, vals.length);
    assertEquals(3, vals[0]);
    assertEquals(1, vals[1]);
    assertEquals("Derived annotation AC missing required fields in VCF header (FORMAT fields: GT)", ann.checkHeader(null));
  }

  public void testNoAlt() {
    final AlleleCountInGenotypesAnnotation ann = new AlleleCountInGenotypesAnnotation();
    final VcfRecord record = new VcfRecord();
    record.setNumberOfSamples(4);
    record.setRefCall("G");
    record.addFormatAndSample("GT", "0/0");
    record.addFormatAndSample("GT", ".");
    assertNull(ann.getValue(record, -1));
  }

}
