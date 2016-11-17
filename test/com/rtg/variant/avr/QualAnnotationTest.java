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

package com.rtg.variant.avr;

import com.rtg.vcf.VcfRecord;

import junit.framework.TestCase;

/**
 */
public class QualAnnotationTest extends TestCase {

  public void test() {
    final QualAnnotation qual = new QualAnnotation();
    assertEquals("QUAL", qual.getName());
    assertEquals(AnnotationDataType.DOUBLE, qual.getType());
    final VcfRecord rec = new VcfRecord("ref", 2, "A");
    assertNull(qual.getValue(rec, -1));
    rec.setQuality(".");
    assertNull(qual.getValue(rec, -1));
    rec.setQuality("10.5");
    final Object val = qual.getValue(rec, 123);
    assertTrue(val instanceof Double);
    assertEquals(10.5, val);
  }
}
