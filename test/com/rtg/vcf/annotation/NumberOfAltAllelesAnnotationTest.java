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
public class NumberOfAltAllelesAnnotationTest extends TestCase {

  public void testName() {
    final NumberOfAltAllelesAnnotation ann = new NumberOfAltAllelesAnnotation();
    assertEquals("NAA", ann.getName());
    assertEquals("Number of alternative alleles", ann.getDescription());
    assertEquals(AnnotationDataType.INTEGER, ann.getType());
  }

  public void test() {
    final NumberOfAltAllelesAnnotation ann = new NumberOfAltAllelesAnnotation();
    VcfRecord rec = new VcfRecord();
    rec.setRefCall("A");
    rec.addAltCall(".");
    assertEquals(0, ann.getValue(rec, 0).intValue());
    rec.setRefCall("A");
    rec.addAltCall("T");
    assertEquals(1, ann.getValue(rec, -2).intValue());
    rec = new VcfRecord();
    rec.setRefCall("A");
    rec.addAltCall("AAAA");
    assertEquals(1, ann.getValue(rec, 1).intValue());
    rec = new VcfRecord();
    rec.setRefCall("AA");
    rec.addAltCall("AAAA");
    rec.addAltCall("AAAAA");
    rec.addAltCall("A");
    assertEquals(3, ann.getValue(rec, 1).intValue());
    rec = new VcfRecord();
    rec.setRefCall(".");
    rec.addAltCall("AAAA");
    rec.addAltCall("AAAAA");
    rec.addAltCall("A");
    assertEquals(3, ann.getValue(rec, 0).intValue());
    assertNull(ann.checkHeader(null));
  }

}
