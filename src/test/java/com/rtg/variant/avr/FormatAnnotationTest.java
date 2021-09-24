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

import com.rtg.vcf.VcfReaderTest;
import com.rtg.vcf.VcfRecord;

import junit.framework.TestCase;

/**
 */
public class FormatAnnotationTest extends TestCase {

  public void testConstructor() {
    try {
      new FormatAnnotation(null, AnnotationDataType.INTEGER);
      fail("accepted null name");
    } catch (NullPointerException npe) {
      // expected
    }

    try {
      new FormatAnnotation("I", null);
      fail("accepted null type");
    } catch (NullPointerException npe) {
      // expected
    }

    final FormatAnnotation fa = new FormatAnnotation("I", AnnotationDataType.INTEGER);
    assertEquals("FORMAT-I", fa.getName());
    assertEquals(AnnotationDataType.INTEGER, fa.getType());
  }

  public void testGetValue() {
    final VcfRecord record = VcfReaderTest.vcfLineToRecord("chr5\t12041\trs55926606\tA\tT\t100\tPASS\tXRX\tGT:I:M:S:F:V\t0|1:12:.:BLAH:1.23:V,V\t.\t1|1:99:.:FART:99.99:.");
    FormatAnnotation fa = new FormatAnnotation("I", AnnotationDataType.INTEGER);
    assertEquals("FORMAT-I", fa.getName());
    assertEquals(AnnotationDataType.INTEGER, fa.getType());
    Object value = fa.getValue(record, 0);
    assertTrue(value instanceof Integer);
    assertEquals(12, value);

    fa = new FormatAnnotation("S", AnnotationDataType.STRING);
    assertEquals("FORMAT-S", fa.getName());
    assertEquals(AnnotationDataType.STRING, fa.getType());
    value = fa.getValue(record, 0);
    assertTrue(value instanceof String);
    assertEquals("BLAH", value);

    fa = new FormatAnnotation("F", AnnotationDataType.DOUBLE);
    assertEquals("FORMAT-F", fa.getName());
    assertEquals(AnnotationDataType.DOUBLE, fa.getType());
    value = fa.getValue(record, 0);
    assertTrue(value instanceof Double);
    assertEquals(1.23, value);

    fa = new FormatAnnotation("M", AnnotationDataType.DOUBLE);
    assertEquals("FORMAT-M", fa.getName());
    assertEquals(AnnotationDataType.DOUBLE, fa.getType());
    value = fa.getValue(record, 0);
    assertNull(value);

    fa = new FormatAnnotation("S", AnnotationDataType.DOUBLE);
    assertEquals("FORMAT-S", fa.getName());
    assertEquals(AnnotationDataType.DOUBLE, fa.getType());
    try {
      fa.getValue(record, 0);
      fail("parsed non number string");
    } catch (NumberFormatException nfe) {
      // expected
    }

    fa = new FormatAnnotation("V", AnnotationDataType.STRING);
    assertEquals("FORMAT-V", fa.getName());
    assertEquals(AnnotationDataType.STRING, fa.getType());
    try {
      fa.getValue(record, 0);
      fail("allowed multiple value field");
    } catch (IllegalArgumentException nfe) {
      // expected
    }

    // sample 1
    fa = new FormatAnnotation("I", AnnotationDataType.INTEGER);
    value = fa.getValue(record, 1);
    assertNull(value);

    fa = new FormatAnnotation("S", AnnotationDataType.STRING);
    value = fa.getValue(record, 1);
    assertNull(value);

    fa = new FormatAnnotation("F", AnnotationDataType.DOUBLE);
    value = fa.getValue(record, 1);
    assertNull(value);

    fa = new FormatAnnotation("M", AnnotationDataType.DOUBLE);
    value = fa.getValue(record, 1);
    assertNull(value);

    // sample 2
    fa = new FormatAnnotation("I", AnnotationDataType.INTEGER);
    value = fa.getValue(record, 2);
    assertTrue(value instanceof Integer);
    assertEquals(99, value);

    fa = new FormatAnnotation("S", AnnotationDataType.STRING);
    value = fa.getValue(record, 2);
    assertTrue(value instanceof String);
    assertEquals("FART", value);

    fa = new FormatAnnotation("F", AnnotationDataType.DOUBLE);
    value = fa.getValue(record, 2);
    assertTrue(value instanceof Double);
    assertEquals(99.99, value);

    fa = new FormatAnnotation("M", AnnotationDataType.INTEGER);
    value = fa.getValue(record, 2);
    assertNull(value);

  }
}
