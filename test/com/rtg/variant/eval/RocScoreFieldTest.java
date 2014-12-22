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
package com.rtg.variant.eval;

import com.rtg.util.TestUtils;
import com.rtg.vcf.VcfRecord;

import junit.framework.TestCase;

/**
 */
public class RocScoreFieldTest extends TestCase {

  public void test() {
    TestUtils.testEnum(RocScoreField.class, "[QUAL, INFO, FORMAT, DERIVED]");
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
    .addInfo("TEST", "45", "46", "47", "48")
    .setNumberOfSamples(2)
    .addFormatAndSample("GT", "0/0")
    .addFormatAndSample("GT", "0/1")
    .addFormatAndSample("GQ", "100")
    .addFormatAndSample("GQ", "95");
    assertEquals(100.0, RocScoreField.FORMAT.getExtractor("GQ", RocSortOrder.DESCENDING).getSortValue(rec, 0));
    assertEquals(95.0, RocScoreField.FORMAT.getExtractor("GQ", RocSortOrder.ASCENDING).getSortValue(rec, 1));
    assertEquals(12.8, RocScoreField.QUAL.getExtractor("GQ", RocSortOrder.DESCENDING).getSortValue(rec, 1));
    assertEquals(23.0, RocScoreField.INFO.getExtractor("DP", RocSortOrder.DESCENDING).getSortValue(rec, 0));
    assertEquals(0.271, RocScoreField.DERIVED.getExtractor("EP", RocSortOrder.DESCENDING).getSortValue(rec, 1), 0.01);
    assertEquals(12.8 / 23.0, RocScoreField.DERIVED.getExtractor("QD", RocSortOrder.DESCENDING).getSortValue(rec, 0), 0.01);
    assertEquals(0.0, RocSortValueExtractor.NULL_EXTRACTOR.getSortValue(rec, 0));
  }
}
