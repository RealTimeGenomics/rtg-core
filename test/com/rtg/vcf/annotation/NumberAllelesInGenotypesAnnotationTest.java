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
public class NumberAllelesInGenotypesAnnotationTest extends TestCase {

  public void testName() {
    final NumberAllelesInGenotypesAnnotation ann = new NumberAllelesInGenotypesAnnotation();
    assertEquals("AN", ann.getName());
    assertEquals("Total number of alleles in called genotypes", ann.getDescription());
    assertEquals(AnnotationDataType.INTEGER, ann.getType());
  }

  public void test() {
    final NumberAllelesInGenotypesAnnotation ann = new NumberAllelesInGenotypesAnnotation();
    VcfRecord record = new VcfRecord();
    record.setNumberOfSamples(3);
    record.setRefCall("G");
    record.addAltCall("A");
    record.addAltCall("C");
    assertEquals(null, ann.getValue(record, -1));

    record.addFormatAndSample("GT", "1/0");
    record.addFormatAndSample("GT", ".");
    record.addFormatAndSample("GT", "2");
    assertEquals(3, ((Integer) ann.getValue(record, -1)).intValue());
    assertEquals("Derived annotation AN missing required fields in VCF header (FORMAT fields: GT)", ann.checkHeader(null));
  }

}
