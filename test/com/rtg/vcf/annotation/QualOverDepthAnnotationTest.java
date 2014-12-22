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
public class QualOverDepthAnnotationTest extends TestCase {
  public void testName() {
    final QualOverDepthAnnotation ann = new QualOverDepthAnnotation();
    assertEquals("QD", ann.getName());
    assertEquals("QUAL / DP", ann.getDescription());
    assertEquals(AnnotationDataType.DOUBLE, ann.getType());
  }

  public void test() {
    final QualOverDepthAnnotation ann = new QualOverDepthAnnotation();
    VcfRecord rec = new VcfRecord();
    rec.setQuality("987.6");
    rec.addInfo("DP", "123");
    assertEquals(8.029, ann.getValue(rec, -1), 0.001);

    rec = new VcfRecord();
    rec.setQuality("987.6");
    rec.addInfo("DP", "0");
    assertEquals(Double.POSITIVE_INFINITY, ann.getValue(rec, -1), 0.001);

    rec = new VcfRecord();
    rec.addInfo("DP", "123");
    assertNull(ann.getValue(rec, 23));

    rec = new VcfRecord();
    rec.setQuality("987.6");
    assertNull(ann.getValue(rec, 0));

    rec = new VcfRecord();
    rec.setQuality("987.6");
    rec.addFormatAndSample("DP", "123");
    assertEquals(8.029, ann.getValue(rec, 0), 0.001);

    rec = new VcfRecord();
    rec.setQuality("987.6");
    rec.addFormatAndSample("DP", "0");
    assertEquals(Double.POSITIVE_INFINITY, ann.getValue(rec, -1), 0.001);

    rec = new VcfRecord();
    rec.setQuality("987.6");
    rec.setNumberOfSamples(2);
    rec.addFormatAndSample("DP", "63");
    rec.addFormatAndSample("DP", "20");
    assertEquals(11.899, ann.getValue(rec, 0), 0.001);

    rec = new VcfRecord();
    rec.setQuality("987.6");
    rec.setNumberOfSamples(3);
    rec.addFormatAndSample("DP", "63");
    rec.addFormatAndSample("DP", "20");
    assertEquals(11.899, ann.getValue(rec, 0), 0.001);

    rec = new VcfRecord();
    rec.setQuality("987.6");
    rec.setNumberOfSamples(3);
    rec.addFormatAndSample("DP", "63");
    rec.addFormatAndSample("DP", ".");
    rec.addFormatAndSample("DP", "20");
    assertEquals(11.899, ann.getValue(rec, 0), 0.001);

    rec = new VcfRecord();
    rec.setQuality("987.6");
    rec.setNumberOfSamples(3);
    rec.addFormatAndSample("DP", "63");
    rec.addFormatAndSample("DP", "0");
    rec.addFormatAndSample("DP", "20");
    assertEquals(11.899, ann.getValue(rec, 0), 0.001);

    rec = new VcfRecord();
    rec.setQuality("987.6");
    rec.addInfo("DP", "100");
    rec.setNumberOfSamples(3);
    rec.addFormatAndSample("DP", "63");
    rec.addFormatAndSample("DP", "0");
    rec.addFormatAndSample("DP", "20");
    assertEquals(9.876, ann.getValue(rec, 0), 0.001);

    assertEquals("Derived annotation QD missing required fields in VCF header (INFO fields: DP)", ann.checkHeader(null));
    final VcfHeader header = new VcfHeader();
    header.addFormatField("DP", MetaType.INTEGER, new VcfNumber("1"), "Depth");
    assertNull(ann.checkHeader(header));
  }

}
