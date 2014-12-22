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

package com.rtg.vcf;

import com.rtg.vcf.annotation.AbstractDerivedAnnotation;
import com.rtg.vcf.annotation.AnnotationDataType;
import com.rtg.vcf.header.MetaType;
import com.rtg.vcf.header.VcfHeader;

import junit.framework.TestCase;

/**
 */
public class VcfInfoDoubleAnnotatorTest extends TestCase {

  private static final class DummyDerivedAnnotation extends AbstractDerivedAnnotation {

    public DummyDerivedAnnotation() {
      super("DUMMY", "DUMMY-DESC", AnnotationDataType.DOUBLE);
    }

    @Override
    public Object getValue(VcfRecord record, int sampleNumber) {
      assertEquals(-1, sampleNumber);
      return 0.1234;
    }

    @Override
    public String checkHeader(VcfHeader header) {
      return null;
    }
  }

  public void test() {
    final VcfInfoDoubleAnnotator ann = new VcfInfoDoubleAnnotator(new DummyDerivedAnnotation());
    final VcfHeader header = new VcfHeader();
    ann.updateHeader(header);
    ann.updateHeader(header); // Test that doing a second time doesn't break it / add extra
    assertEquals(1, header.getInfoLines().size());
    assertEquals("DUMMY", header.getInfoLines().get(0).getId());
    assertEquals("DUMMY-DESC", header.getInfoLines().get(0).getDescription());
    assertEquals(MetaType.FLOAT, header.getInfoLines().get(0).getType());
    assertEquals(1, header.getInfoLines().get(0).getNumber().getNumber());
    final VcfRecord rec = new VcfRecord();
    ann.annotate(rec);
    ann.annotate(rec); // Test that doing a second time doesn't break it / add extra
    assertEquals(1, rec.getInfo().get("DUMMY").size());
    assertEquals("0.123", rec.getInfo().get("DUMMY").get(0));
    assertEquals("null\t0\t.\tnull\t.\t.\t.\tDUMMY=0.123", rec.toString());
  }
}
