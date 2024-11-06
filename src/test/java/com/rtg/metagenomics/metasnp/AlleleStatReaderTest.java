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

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.util.Arrays;

import com.rtg.util.StringUtils;

import junit.framework.TestCase;

/**
 */
public class AlleleStatReaderTest extends TestCase {

  private static final String SIMPLE = "sequence\tposition\treference\ta\tc\tg\tt\taRatio\tcRatio\tgRatio\ttRatio" + StringUtils.LS
    + "seq\t1\tC\t0,0,0,0\t88,122,104,134\t0,0,0,0\t0,0,0,0\t0.000000,0.000000,0.000000,0.000000\t-1,-1,-1,-1\t0.000000,0.000000,0.000000,0.000000\t0.000000,0.000000,0.000000,0.000000" + StringUtils.LS
    + "seq\t2\tG\t0,0,0,0\t0,0,0,0\t93,127,109,137\t0,0,0,0\t0.000000,0.000000,0.000000,0.000000\t0.000000,0.000000,0.000000,0.000000\t-1,-1,-1,-1\t0.000000,0.000000,0.000000,0.000000" + StringUtils.LS
    + "seq\t3\tA\t99,130,112,148\t0,0,0,0\t0,0,0,0\t0,0,0,0\t-1,-1,-1,-1\t0.000000,0.000000,0.000000,0.000000\t0.000000,0.000000,0.000000,0.000000\t0.000000,0.000000,0.000000,0.000000" + StringUtils.LS
    + "seq2\t46\tT\t0,0,0,0\t0,0,0,0\t0,0,0,0\t102,135,117,155" + StringUtils.LS;

  public void testSimple() throws IOException {
    try (AlleleStatReader reader = new AlleleStatReader(new ByteArrayInputStream(SIMPLE.getBytes()))) {
      final MetaSnpLine line = reader.nextLine();
      assertEquals("seq", line.getSequence());
      assertEquals(0, line.getPosition());
      assertEquals(1, line.getReferenceIndex());
      assertTrue(Arrays.equals(new double[]{0, 0, 0, 0}, line.mCounts[0]));
      assertTrue(Arrays.equals(new double[]{88, 122, 104, 134}, line.mCounts[1]));
      assertNotNull(reader.nextLine());
      assertNotNull(reader.nextLine());

      final MetaSnpLine line2 = reader.nextLine();
      assertEquals("seq2", line2.getSequence());
      assertEquals(45, line2.getPosition());
      assertEquals(3, line2.getReferenceIndex());
      assertEquals("T", line2.getReferenceAllele());
      assertTrue(Arrays.equals(new double[]{0, 0, 0, 0}, line2.mCounts[0]));
      assertTrue(Arrays.equals(new double[]{102, 135, 117, 155}, line2.mCounts[3]));
      assertNull(reader.nextLine());
    }
  }
}
