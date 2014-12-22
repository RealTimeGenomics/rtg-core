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
import com.rtg.vcf.header.VcfHeader;

import junit.framework.TestCase;

/**
 */
public class DummyDerivedFormatAnnotationTest extends TestCase {

  private static final class DummyDerivedFormatAnnotation extends AbstractDerivedFormatAnnotation {
    public DummyDerivedFormatAnnotation() {
      super("DUMMY", "DUMMY-DESC", AnnotationDataType.STRING);
    }
    @Override
    public Object getValue(VcfRecord record, int sampleNumber) {
      return sampleNumber;
    }
    @Override
    public String checkHeader(VcfHeader header) {
      return null;
    }
  }

  public void testDummy() {
    final AbstractDerivedFormatAnnotation anno = new DummyDerivedFormatAnnotation();
    assertEquals("DUMMY", anno.getName());
    assertEquals("DUMMY-DESC", anno.getDescription());
    assertEquals(1, anno.getValue(null, 1));
  }
}
