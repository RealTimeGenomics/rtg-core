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

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;

import com.rtg.util.TestUtils;
import com.rtg.vcf.annotation.DerivedAnnotations;

import junit.framework.TestCase;

/**
 */
public class AnnotationLoaderTest extends TestCase {

  public void testEnum() {
    TestUtils.testEnum(AnnotationLoader.AnnotationType.class, "[QUAL, FORMAT, INFO, DERIVED]");
    assertEquals(0, AnnotationLoader.AnnotationType.QUAL.ordinal());
    assertEquals(1, AnnotationLoader.AnnotationType.FORMAT.ordinal());
    assertEquals(2, AnnotationLoader.AnnotationType.INFO.ordinal());
    assertEquals(3, AnnotationLoader.AnnotationType.DERIVED.ordinal());

    // new annotations need to be added at end of enum
    assertEquals(4, AnnotationLoader.AnnotationType.values().length);
  }

  public void testLoad() throws IOException {
    try (final ByteArrayOutputStream baos = new ByteArrayOutputStream()) {
      try (final DataOutputStream dos = new DataOutputStream(baos)) {
        final Annotation[] annos = {
            new QualAnnotation(), new InfoAnnotation("info", AnnotationDataType.STRING), new FormatAnnotation("format", AnnotationDataType.INTEGER), new DerivedAnnotation(DerivedAnnotations.IC.getAnnotation())
        };
        for (Annotation anno : annos) {
          anno.save(dos);
        }
      }
      try (final DataInputStream dis = new DataInputStream(new ByteArrayInputStream(baos.toByteArray()))) {
        Annotation anno = AnnotationLoader.load(dis);
        assertTrue(anno instanceof QualAnnotation);
        anno = AnnotationLoader.load(dis);
        assertTrue(anno instanceof InfoAnnotation);
        assertEquals("INFO-info", anno.getName());
        assertEquals(AnnotationDataType.STRING, anno.getType());
        anno = AnnotationLoader.load(dis);
        assertTrue(anno instanceof FormatAnnotation);
        assertEquals(AnnotationDataType.INTEGER, anno.getType());
        assertEquals("FORMAT-format", anno.getName());
        anno = AnnotationLoader.load(dis);
        assertTrue(anno instanceof DerivedAnnotation);
        assertEquals(AnnotationDataType.DOUBLE, anno.getType());
        assertEquals("DERIVED-IC", anno.getName());
      }
    }
  }

}
