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

package com.rtg.metagenomics.metasnp;

import java.io.IOException;

import com.rtg.vcf.VcfRecord;

import junit.framework.TestCase;

/**
 * Tests the corresponding class.
 */
public class MetaSnpLineTest extends TestCase {

  public void testAlleleStats() throws IOException {
    final MetaSnpLine line = MetaSnpLine.create("chr1\t3\tc\t1.9\t6.0\t0\t0", 0);
    assertEquals("chr1", line.getSequence());
    assertEquals(2, line.getPosition()); // 0-based
    assertEquals("c", line.getReferenceAllele());
    assertEquals(4, line.mCounts.length);
    assertEquals(1.9, line.mCounts[0][0], 1e-9);
    assertEquals(6.0, line.mCounts[1][0], 1e-9);
    assertEquals(0, line.mCounts[2][0], 1e-9);
    assertEquals(0, line.mCounts[3][0], 1e-9);
  }

  public void testVcfRecord() {
    final VcfRecord rec = new VcfRecord("chr1", 2, "c"); // position is 0-based in constructor
    rec.setId(".")
      .setQuality("12.8")
      .addAltCall("a")
      .setNumberOfSamples(1)
      .addFormatAndSample("AD", "6.0,1.9");
    final MetaSnpLine line = MetaSnpLine.create(rec);
    assertEquals("chr1", line.getSequence());
    assertEquals(2, line.getPosition()); // 0-based
    assertEquals("c", line.getReferenceAllele());
    assertEquals(2, line.mCounts.length);
    assertEquals(6.0, line.mCounts[0][0], 1e-9); // ref allele, c in this case
    assertEquals(1.9, line.mCounts[1][0], 1e-9); // first alt allele, a in this case
  }
}
