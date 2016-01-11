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

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.util.Arrays;

import com.rtg.util.StringUtils;

import junit.framework.TestCase;

/**
 */
public class AlleleStatReaderTest extends TestCase {
static final String SIMPLE = "sequence\tposition\treference\ta\tc\tg\tt\taRatio\tcRatio\tgRatio\ttRatio" + StringUtils.LS
                             + "seq\t1\tC\t0,0,0,0\t88,122,104,134\t0,0,0,0\t0,0,0,0\t0.000000,0.000000,0.000000,0.000000\t-1,-1,-1,-1\t0.000000,0.000000,0.000000,0.000000\t0.000000,0.000000,0.000000,0.000000" + StringUtils.LS
                             + "seq\t2\tG\t0,0,0,0\t0,0,0,0\t93,127,109,137\t0,0,0,0\t0.000000,0.000000,0.000000,0.000000\t0.000000,0.000000,0.000000,0.000000\t-1,-1,-1,-1\t0.000000,0.000000,0.000000,0.000000" + StringUtils.LS
                             + "seq\t3\tA\t99,130,112,148\t0,0,0,0\t0,0,0,0\t0,0,0,0\t-1,-1,-1,-1\t0.000000,0.000000,0.000000,0.000000\t0.000000,0.000000,0.000000,0.000000\t0.000000,0.000000,0.000000,0.000000" + StringUtils.LS
                             + "seq2\t46\tT\t0,0,0,0\t0,0,0,0\t0,0,0,0\t102,135,117,155" + StringUtils.LS;

  public void testSimple() throws IOException {
    try (AlleleStatReader reader = new AlleleStatReader(new ByteArrayInputStream(SIMPLE.getBytes()))) {
      final MetaSnpLine line = reader.nextLine();
      assertEquals("seq", line.getSequence());
      assertEquals(0, line.getPosition());
      assertEquals(1, line.mReference);
      assertTrue(Arrays.equals(new double[] {0, 0, 0, 0}, line.mCounts[0]));
      assertTrue(Arrays.equals(new double[] {88, 122, 104, 134}, line.mCounts[1]));
      assertNotNull(reader.nextLine());
      assertNotNull(reader.nextLine());

      final MetaSnpLine line2 = reader.nextLine();
      assertEquals("seq2", line2.getSequence());
      assertEquals(45, line2.getPosition());
      assertEquals(3, line2.mReference);
      assertTrue(Arrays.equals(new double[] {0, 0, 0, 0}, line2.mCounts[0]));
      assertTrue(Arrays.equals(new double[] {102, 135, 117, 155}, line2.mCounts[3]));
      assertNull(reader.nextLine());
    }

  }
}
