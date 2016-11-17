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

import com.rtg.vcf.VcfReader;
import com.rtg.vcf.VcfRecord;

import junit.framework.TestCase;

/**
 */
public class InfoAnnotationTest extends TestCase {

  public void testConstructor() {
    try {
      new InfoAnnotation(null, AnnotationDataType.INTEGER);
      fail("accepted null name");
    } catch (NullPointerException npe) {
      // expected
    }

    final InfoAnnotation fa = new InfoAnnotation("I", AnnotationDataType.INTEGER);
    assertEquals("INFO-I", fa.getName());
    assertEquals(AnnotationDataType.INTEGER, fa.getType());
  }

  public void testGetValue() {
    final VcfRecord record = VcfReader.vcfLineToRecord("chr5\t12041\trs55926606\tA\tT\t100\tPASS\tXRX;IS=BOO;IF=98.7;II=42;CC=a,b;DD=D1\tGT:I:M:S:F\t0|1:12:.:BLAH:1.23\t.\t1|1:99:.:FART:99.99");
    InfoAnnotation fa = new InfoAnnotation("XRX", AnnotationDataType.BOOLEAN);
    assertEquals("INFO-XRX", fa.getName());
    assertEquals(AnnotationDataType.BOOLEAN, fa.getType());
    Object value = fa.getValue(record, -1);
    assertTrue(value instanceof Boolean);
    assertEquals(Boolean.TRUE, value);

    fa = new InfoAnnotation("YRY", AnnotationDataType.BOOLEAN);
    assertEquals("INFO-YRY", fa.getName());
    assertEquals(AnnotationDataType.BOOLEAN, fa.getType());
    value = fa.getValue(record, 0);
    assertTrue(value instanceof Boolean);
    assertEquals(Boolean.FALSE, value);

    fa = new InfoAnnotation("IS", AnnotationDataType.STRING);
    assertEquals("INFO-IS", fa.getName());
    assertEquals(AnnotationDataType.STRING, fa.getType());
    value = fa.getValue(record, 123);
    assertTrue(value instanceof String);
    assertEquals("BOO", value);

    fa = new InfoAnnotation("IF", AnnotationDataType.DOUBLE);
    assertEquals("INFO-IF", fa.getName());
    assertEquals(AnnotationDataType.DOUBLE, fa.getType());
    value = fa.getValue(record, 0);
    assertTrue(value instanceof Double);
    assertEquals(98.7, value);

    fa = new InfoAnnotation("II", AnnotationDataType.INTEGER);
    assertEquals("INFO-II", fa.getName());
    assertEquals(AnnotationDataType.INTEGER, fa.getType());
    value = fa.getValue(record, -1);
    assertTrue(value instanceof Integer);
    assertEquals(42, value);

    fa = new InfoAnnotation("IS", AnnotationDataType.DOUBLE);
    assertEquals("INFO-IS", fa.getName());
    assertEquals(AnnotationDataType.DOUBLE, fa.getType());
    try {
      fa.getValue(record, 0);
      fail("parsed non number string");
    } catch (NumberFormatException nfe) {
      // expected
    }

    fa = new InfoAnnotation("XRX", AnnotationDataType.INTEGER);
    assertEquals("INFO-XRX", fa.getName());
    assertEquals(AnnotationDataType.INTEGER, fa.getType());
    try {
      fa.getValue(record, 1);
      fail("allowed no value");
    } catch (IllegalArgumentException iae) {
      // expected
    }

    fa = new InfoAnnotation("CC", AnnotationDataType.STRING);
    assertEquals("INFO-CC", fa.getName());
    assertEquals(AnnotationDataType.STRING, fa.getType());
    try {
      fa.getValue(record, 1);
      fail("allowed more than one value");
    } catch (IllegalArgumentException iae) {
      // expected
    }
  }

}
