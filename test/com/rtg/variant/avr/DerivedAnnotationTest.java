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
import com.rtg.vcf.annotation.AbstractDerivedInfoAnnotation;
import com.rtg.vcf.header.InfoField;
import com.rtg.vcf.header.MetaType;
import com.rtg.vcf.header.VcfHeader;
import com.rtg.vcf.header.VcfNumber;

import junit.framework.TestCase;

/**
 */
public class DerivedAnnotationTest extends TestCase {

  public void testConstructor() {
    final DerivedAnnotation ann = new DerivedAnnotation("IC");
    assertEquals("DERIVED-IC", ann.getName());
    assertEquals(AnnotationDataType.DOUBLE, ann.getType());
  }

  private static final class DummyVcfAnnotation extends AbstractDerivedInfoAnnotation {
    DummyVcfAnnotation() {
      super(new InfoField("DUMMY", MetaType.FLAG, VcfNumber.ONE, "DUMMY-DESC"), null);
    }
    @Override
    public Object getValue(VcfRecord record, int sampleNumber) {
      assertNull(record);
      return true;
    }
    @Override
    public String checkHeader(VcfHeader header) {
      assertNull(header);
      return "A String";
    }
  }

  public void testWrapping() {
    final DerivedAnnotation ann = new DerivedAnnotation(new DummyVcfAnnotation());
    assertEquals("DERIVED-DUMMY", ann.getName());
    assertEquals(AnnotationDataType.BOOLEAN, ann.getType());
    assertEquals(true, ann.getValue(null, 0));
    assertEquals("A String", ann.checkHeader(null));
  }
}
