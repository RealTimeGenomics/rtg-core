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

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;

import com.rtg.util.TestUtils;
import com.rtg.vcf.annotation.AnnotationDataType;
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
