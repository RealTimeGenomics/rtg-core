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

import static com.rtg.util.StringUtils.TAB;

import java.io.ByteArrayOutputStream;
import java.io.IOException;

import com.rtg.vcf.header.VcfHeader;

import junit.framework.TestCase;

/**
 */
public class VcfWriterTest extends TestCase {

  private static final String LS = "\n"; // Not platform specific

  public VcfWriterTest(String name) {
    super(name);
  }

  public void test() throws IOException {
    final VcfHeader head = new VcfHeader();
    head.setVersionValue(VcfHeader.VERSION_VALUE);
    head.addMetaInformationLine("##test1212121")
    .addMetaInformationLine("##test12");
    head.addSampleName("sample1")
    .addSampleName("sample2");

    final VcfRecord rec = new VcfRecord();
    rec.setSequence("chr1")
    .setStart(1209)
    .setId(".")
    .setQuality("12.8")
    .setRefCall("a")
    .addAltCall("c")
    .addAltCall("t")
    .addFilter("TEST1")
    .addFilter("TEST2")
    .addInfo("DP", "23")
    .addInfo("TEST", "45,46,47,48")
    .setNumberOfSamples(2)
    .addFormatAndSample("GT", "0/0")
    .addFormatAndSample("GT", "0/1")
    .addFormatAndSample("GQ", "100")
    .addFormatAndSample("GQ", "95")
    ;

    final ByteArrayOutputStream bos = new ByteArrayOutputStream();
    final VcfWriter w = new VcfWriter(head, bos);
    w.write(rec);
    w.write(rec);

    final String line = ""
      + "chr1" + TAB
      + "1210" + TAB
      + "." + TAB
      + "a" + TAB
      + "c,t" + TAB
      + "12.8" + TAB
      + "TEST1;TEST2" + TAB
      + "DP=23;TEST=45,46,47,48" + TAB
      + "GT:GQ" + TAB
      + "0/0:100" + TAB
      + "0/1:95" + LS;
    final String exp = "##fileformat=VCFv4.1" + LS
      + "##test1212121" + LS
      + "##test12" + LS
      + "#CHROM" + TAB + "POS" + TAB + "ID" + TAB + "REF" + TAB + "ALT" + TAB + "QUAL" + TAB + "FILTER" + TAB + "INFO" + TAB + "FORMAT" + TAB + "sample1" + TAB + "sample2" + LS
      + line
      + line;

    assertEquals(exp, bos.toString());
    assertEquals(w.getHeader(), head);
  }

  public void testErrors() {
    try {
      new VcfWriter(null, new ByteArrayOutputStream());
    } catch (NullPointerException ex) {
      assertEquals("header cannot be null", ex.getMessage());
    }

    try {
      new VcfWriter(new VcfHeader(), null);
    } catch (NullPointerException ex) {
      assertEquals("output stream cannot be null", ex.getMessage());
    }
  }
}
