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
public class LongestAlleleAnnotationTest extends TestCase {

  public void testName() {
    final LongestAlleleAnnotation lalAnn = new LongestAlleleAnnotation();
    assertEquals("LAL", lalAnn.getName());
    assertEquals("Length of longest allele", lalAnn.getDescription());
    assertEquals(AnnotationDataType.INTEGER, lalAnn.getType());
  }

  public void test() {
    final LongestAlleleAnnotation lalAnn = new LongestAlleleAnnotation();
    VcfRecord rec = new VcfRecord();
    rec.setRefCall("AAA");
    rec.addAltCall("A");
    assertEquals(3, ((Integer) lalAnn.getValue(rec, 0)).intValue());
    rec = new VcfRecord();
    rec.setRefCall("A");
    rec.addAltCall("AAAA");
    assertEquals(4, ((Integer) lalAnn.getValue(rec, -1)).intValue());
    rec = new VcfRecord();
    rec.setRefCall("AA");
    rec.addAltCall("AAAA");
    rec.addAltCall("AAAAA");
    rec.addAltCall("A");
    assertEquals(5, ((Integer) lalAnn.getValue(rec, 34)).intValue());
    assertNull(lalAnn.checkHeader(null));
  }
}
