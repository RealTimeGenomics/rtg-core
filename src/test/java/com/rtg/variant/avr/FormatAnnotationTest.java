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
