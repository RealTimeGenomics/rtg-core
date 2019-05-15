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

package com.rtg.visualization;

import static com.rtg.util.StringUtils.TAB;

import com.rtg.util.MathUtils;
import com.rtg.util.PosteriorUtils;
import com.rtg.vcf.VcfReaderTest;
import com.rtg.vcf.VcfRecord;

import junit.framework.TestCase;

/**
 */
public class AviewVariantTest extends TestCase {

  public void test() {
    final String vcf = "someKindOfName" + TAB + "23" + TAB + "." + TAB + "A" + TAB + "T" + TAB + "12.8" + TAB + "PASS" + TAB + "." + TAB + "GT" + TAB + "1/1";
    final VcfRecord rec = VcfReaderTest.vcfLineToRecord(vcf);

    final AviewVariant v = new AviewVariant(rec, 0);

    assertEquals(23, v.getPosition());
    assertEquals(1, v.referenceLength());
  }

  private AviewVariant getAviewVariant(String var) {
    return new AviewVariant(VcfReaderTest.vcfLineToRecord(var.replaceAll(" ", TAB)), 0);
  }

  private static final String SNP_LINE = "someKindOfName" + TAB + "23" + TAB + "." + TAB + "A" + TAB + "T" + TAB + "12.8" + TAB + "PASS" + TAB + "." + TAB + "GT:DP:RE:GQ" + TAB + "1/1:4:0.02:" + PosteriorUtils.phredIfy(12.8 * MathUtils.LOG_10);
  public void testSnpConstruction() {
    final AviewVariant variant = getAviewVariant(SNP_LINE);
    assertEquals(23, variant.getPosition());
    assertEquals(1, variant.nt(true).length());
    assertEquals("T", variant.nt(true));
    assertNull(variant.nt(false));
  }

  private static final String SNP_LINE2 = "someKindOfName" + TAB + "23" + TAB + "." + TAB + "A" + TAB + "T,C" + TAB + "12.8" + TAB + "PASS" + TAB + "." + TAB + "GT:DP:RE:GQ" + TAB + "1/2:4:0.02:" + PosteriorUtils.phredIfy(12.8 * MathUtils.LOG_10);
  public void testHeterozygousSnpConstruction() {
    final AviewVariant variant = getAviewVariant(SNP_LINE2);
    assertEquals(23, variant.getPosition());
    assertEquals(1, variant.nt(true).length());
    assertEquals("T", variant.nt(true));
    assertEquals(1, variant.nt(false).length());
    assertEquals("C", variant.nt(false));
  }

  private static final String INSERT_LINE = "someKindOfName" + TAB + "22" + TAB + "." + TAB + "A" + TAB + "AACT" + TAB + "12.8" + TAB + "PASS" + TAB + "." + TAB + "GT:DP:RE:GQ" + TAB + "1/1:4:0.02:" + PosteriorUtils.phredIfy(12.8 * MathUtils.LOG_10);
  public void testInsertConstruction() {
    final AviewVariant variant = getAviewVariant(INSERT_LINE);
    assertEquals(23, variant.getPosition());
  }


  private static final String DELETION_LINE = "someKindOfName" + TAB + "22" + TAB + "." + TAB + "GAGT" + TAB + "G" + TAB + "12.8" + TAB + "PASS" + TAB + "." + TAB + "GT:DP:RE:GQ" + TAB + "1/1:4:0.02:" + PosteriorUtils.phredIfy(12.8 * MathUtils.LOG_10);
  public void testDeletionConstructor() {
    final AviewVariant variant = getAviewVariant(DELETION_LINE);
    assertEquals(23, variant.getPosition());
  }

  private static final String MNP_LINE = "someKindOfName" + TAB + "23" + TAB + "." + TAB + "AGT" + TAB + "CTC" + TAB + "12.8" + TAB + "PASS" + TAB + "." + TAB + "GT:DP:RE:GQ" + TAB + "1/1:4:0.02:" + PosteriorUtils.phredIfy(12.8 * MathUtils.LOG_10);
  public void testMnpConstructor() {
    final AviewVariant variant = getAviewVariant(MNP_LINE);
    assertEquals(23, variant.getPosition());

    assertEquals(3, variant.nt(true).length());
    assertEquals("CTC", variant.nt(true));
    assertNull(variant.nt(false));
  }

  private static final String UNCHANGED_LINE = "someKindOfName" + TAB + "23" + TAB + "." + TAB + "A" + TAB + "C" + TAB + "12.8" + TAB + "PASS" + TAB + "." + TAB + "GT:DP:RE:GQ" + TAB + "1/1:4:0.02:" + PosteriorUtils.phredIfy(12.8 * MathUtils.LOG_10);
  public void testUnchangedConstructor() {
    final AviewVariant variant = getAviewVariant(UNCHANGED_LINE);
    assertEquals(23, variant.getPosition());
    assertEquals("C", variant.nt(true));
  }

  private static final String SHORT_LINE = "someKindOfName" + TAB + "23" + TAB + "." + TAB + "A" + TAB + "C" + TAB + "0.0" + TAB + "PASS" + TAB + "." + TAB + "GT" + TAB + "1/1";
  public void testShortConstructor() {
    final AviewVariant variant = getAviewVariant(SHORT_LINE);
    assertEquals(23, variant.getPosition());
    assertEquals("C", variant.nt(true));
  }

  public void testAllFieldsPresent() {
    final String line = ("simulatedSequence1 2180 . C G,T " + PosteriorUtils.phredIfy(45.8) + TAB + "PASS" + TAB + "." + TAB + "GT:DP:RE:GQ" + TAB + "1/2:35:0.697:" + PosteriorUtils.phredIfy(31.0 * MathUtils.LOG_10)).replaceAll(" ", "\t");
    final AviewVariant v = getAviewVariant(line);
    assertEquals(2180, v.getPosition());
    assertEquals(1, v.referenceLength());
    assertEquals("G", v.ntAlleleA());
    assertEquals("T", v.ntAlleleB());
  }

  public void testMissingGT() {
    final String line = "simulatedSequence1 2180 . C T 0.0 PASS . GT .".replaceAll(" ", "\t");
    final AviewVariant v = getAviewVariant(line);
    assertEquals(2180, v.getPosition());
  }

}
