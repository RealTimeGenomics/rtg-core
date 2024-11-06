/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
    assertEquals("C", line.getReferenceAllele());
    assertEquals(4, line.mCounts.length);
    assertEquals(1.9, line.mCounts[0][0], 1e-9);
    assertEquals(6.0, line.mCounts[1][0], 1e-9);
    assertEquals(0, line.mCounts[2][0], 1e-9);
    assertEquals(0, line.mCounts[3][0], 1e-9);
  }

  private void fieldTestVcfRecord(final String ad) {
    final VcfRecord rec = new VcfRecord("chr1", 2, "C"); // position is 0-based in constructor
    rec.setId(".")
      .setQuality("12.8")
      .addAltCall("A")
      .setNumberOfSamples(1)
      .addFormatAndSample(ad, "6.0,1.9");
    final MetaSnpLine line = MetaSnpLine.create(rec);
    assertEquals("chr1", line.getSequence());
    assertEquals(2, line.getPosition()); // 0-based
    assertEquals("C", line.getReferenceAllele());
    assertEquals(2, line.mCounts.length);
    assertEquals(6.0, line.mCounts[0][0], 1e-9); // ref allele, c in this case
    assertEquals(1.9, line.mCounts[1][0], 1e-9); // first alt allele, a in this case
  }

  public void testVcfRecord() {
    fieldTestVcfRecord("ADE");
    fieldTestVcfRecord("AD"); // check that the AD fallback also works
  }
}
