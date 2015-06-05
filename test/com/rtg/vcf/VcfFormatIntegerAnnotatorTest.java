/*
 * Copyright (c) 2015. Real Time Genomics Limited.
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

import com.rtg.vcf.annotation.AbstractDerivedFormatAnnotation;
import com.rtg.vcf.annotation.AnnotationDataType;
import com.rtg.vcf.header.MetaType;
import com.rtg.vcf.header.VcfHeader;

import junit.framework.TestCase;

/**
 * Annotator for integer FORMAT fields.
 */
public class VcfFormatIntegerAnnotatorTest extends TestCase {

  private static final class DummyDerivedAnnotation extends AbstractDerivedFormatAnnotation {
    public DummyDerivedAnnotation() {
      super("DUMMY", "DUMMY-DESC", AnnotationDataType.INTEGER);
    }
    @Override
    public Object getValue(VcfRecord record, int sampleNumber) {
      return sampleNumber * 17;
    }
    @Override
    public String checkHeader(VcfHeader header) {
      return null;
    }
  }

  public void test() {
    final VcfFormatIntegerAnnotator ann = new VcfFormatIntegerAnnotator(new DummyDerivedAnnotation());
    final VcfHeader header = new VcfHeader();
    ann.updateHeader(header);
    ann.updateHeader(header); // Test that doing a second time doesn't break it / add extra
    assertEquals(1, header.getFormatLines().size());
    assertEquals("DUMMY", header.getFormatLines().get(0).getId());
    assertEquals("DUMMY-DESC", header.getFormatLines().get(0).getDescription());
    assertEquals(MetaType.INTEGER, header.getFormatLines().get(0).getType());
    assertEquals(1, header.getFormatLines().get(0).getNumber().getNumber());
    final VcfRecord rec = new VcfRecord();
    rec.setNumberOfSamples(2);
    ann.annotate(rec);
    ann.annotate(rec); // Test that doing a second time doesn't break it / add extra
    assertEquals(2, rec.getFormatAndSample().get("DUMMY").size());
    assertEquals("0", rec.getFormatAndSample().get("DUMMY").get(0));
    assertEquals("17", rec.getFormatAndSample().get("DUMMY").get(1));
    assertEquals("null\t0\t.\tnull\t.\t.\t.\t.\tDUMMY\t0\t17", rec.toString());
  }
}
