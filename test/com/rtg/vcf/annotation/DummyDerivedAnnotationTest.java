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
import com.rtg.vcf.header.MetaType;
import com.rtg.vcf.header.VcfHeader;
import com.rtg.vcf.header.VcfNumber;

import junit.framework.TestCase;

/**
 */
public class DummyDerivedAnnotationTest extends TestCase {

  private static final class DummyDerivedAnnotation extends AbstractDerivedAnnotation {

    public DummyDerivedAnnotation() {
      super("DUMMY", "DUMMY-DESC", AnnotationDataType.DOUBLE);
    }

    @Override
    public Object getValue(VcfRecord record, int sampleNumber) {
      return 0.0;
    }

    @Override
    public String checkHeader(VcfHeader header) {
      return null;
    }
  }

  public void testDerivedAnnotation() {
    final AbstractDerivedAnnotation ann = new DummyDerivedAnnotation();
    assertEquals("DUMMY", ann.getName());
    assertEquals("DUMMY-DESC", ann.getDescription());
    assertEquals(AnnotationDataType.DOUBLE, ann.getType());
  }

  public void testCheckHeader() {
    final AbstractDerivedAnnotation ann = new DummyDerivedAnnotation();
    assertEquals("Derived annotation DUMMY missing required fields in VCF header (INFO fields: II) (FORMAT fields: FF)", ann.checkHeader(null, new String[]{"II"}, new String[] {"FF"}));
    final VcfHeader header = new VcfHeader();
    header.addInfoField("II", MetaType.INTEGER, new VcfNumber("1"), "Info Field");
    header.addFormatField("FF", MetaType.INTEGER, new VcfNumber("1"), "Format Field");
    final String res = ann.checkHeader(header, new String[]{"II"}, new String[] {"FF"});
    assertNull(res, res);
  }
}
