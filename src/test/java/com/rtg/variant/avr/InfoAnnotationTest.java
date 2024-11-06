/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

package com.rtg.variant.avr;

import com.rtg.vcf.VcfReaderTest;
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
    final VcfRecord record = VcfReaderTest.vcfLineToRecord("chr5\t12041\trs55926606\tA\tT\t100\tPASS\tXRX;IS=BOO;IF=98.7;II=42;CC=a,b;DD=D1\tGT:I:M:S:F\t0|1:12:.:BLAH:1.23\t.\t1|1:99:.:FART:99.99");
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
