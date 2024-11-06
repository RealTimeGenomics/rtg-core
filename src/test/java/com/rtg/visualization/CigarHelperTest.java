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
package com.rtg.visualization;

import java.util.Arrays;
import java.util.HashMap;

import htsjdk.samtools.SAMRecord;

import junit.framework.TestCase;

/**
 *
 */
public class CigarHelperTest extends TestCase {


  public void testCigarHelper() {
    final int[] inserts = new int[4];
    CigarHelper.locateInserts("1=1I2X", 0, inserts);
    assertTrue(Arrays.equals(new int[] {0, 1, 0, 0}, inserts));

    Arrays.fill(inserts, 0);
    CigarHelper.locateInserts("4=", 0, inserts);
    assertTrue(Arrays.equals(new int[] {0, 0, 0, 0}, inserts));

    Arrays.fill(inserts, 0);
    CigarHelper.locateInserts("4X", 0, inserts);
    assertTrue(Arrays.equals(new int[] {0, 0, 0, 0}, inserts));

    Arrays.fill(inserts, 0);
    CigarHelper.locateInserts("4M", 0, inserts);
    assertTrue(Arrays.equals(new int[] {0, 0, 0, 0}, inserts));

    Arrays.fill(inserts, 0);
    CigarHelper.locateInserts("4D", 0, inserts);
    assertTrue(Arrays.equals(new int[] {0, 0, 0, 0}, inserts));

    Arrays.fill(inserts, 0);
    CigarHelper.locateInserts("1=3I", 0, inserts);
    assertTrue(Arrays.equals(new int[] {0, 3, 0, 0}, inserts));

    final int[] inserts21 = new int[21];
    CigarHelper.locateInserts("1=13I7X", 0, inserts21);
    //System.err.println(Arrays.toString(inserts20));
    assertTrue(Arrays.equals(new int[] {0, 13, 0, 0, 0, 0, 0, 0, 0, 0,
                                        0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, inserts21));

    Arrays.fill(inserts, 0);
    CigarHelper.locateInserts("1I1N1I", 0, inserts);
    assertTrue(Arrays.equals(new int[] {1, 1, 0, 0}, inserts));

    Arrays.fill(inserts, 0);
    CigarHelper.locateInserts("1S1I1N1I", 0, inserts);
    assertTrue("Exp: [0, 1, 1, 0] Got: " + Arrays.toString(inserts), Arrays.equals(new int[] {1, 1, 0, 0}, inserts));

    Arrays.fill(inserts, 0);
    CigarHelper.locateInserts("1I1N2I", 1, inserts);
    assertTrue(Arrays.equals(new int[] {0, 1, 2, 0}, inserts));

    Arrays.fill(inserts, 0);
    CigarHelper.locateInserts("2=2B2I", 1, inserts);
    assertTrue(Arrays.toString(inserts), Arrays.equals(new int[] {0, 2, 0, 0}, inserts));

    Arrays.fill(inserts, 0);
    CigarHelper.locateInserts("2=2B2=", 1, inserts);
    assertTrue(Arrays.toString(inserts), Arrays.equals(new int[] {0, 0, 0, 0}, inserts));

    Arrays.fill(inserts, 0);
    CigarHelper.locateInserts("1I1X1B2I1=", 1, inserts);
    assertTrue(Arrays.toString(inserts), Arrays.equals(new int[] {0, 2, 0, 0}, inserts));

    Arrays.fill(inserts, 0);
    CigarHelper.locateInserts("2I1X1B1I1=", 1, inserts);
    assertTrue(Arrays.toString(inserts), Arrays.equals(new int[] {0, 2, 0, 0}, inserts));

    Arrays.fill(inserts, 0);
    CigarHelper.locateInserts("2I1X1B2I1=", 1, inserts);
    assertTrue(Arrays.toString(inserts), Arrays.equals(new int[] {0, 2, 0, 0}, inserts));

    final int[] moreInserts = new int[40];
    Arrays.fill(moreInserts, 0);
    CigarHelper.locateInserts("5=2B20=6N6=1I1=1X1=", 0, moreInserts);

    assertEquals(1, moreInserts[35]);

  }

  public void testIlleagalCigar() {
    try {
      CigarHelper.locateInserts("1O", 0, null);
    } catch (final IllegalStateException ex) {
      assertTrue(ex.getMessage().contains("Invalid cigar : 1O"));
    }
  }

  // 1234567890123456789012345678901234567890
  // GACTG
  //     TTGAGAGGGANAAAGTTATG      AACATTTATN

  public void testUnrollLegacyCgCigar() {
    //30873 67  paolo-bac 66  255 4=1X19=6N10=  = 369 303
    SAMRecord r = SamTestHelper.getLegacyCGSamRecord("GACTTTGAGAGGGANAAAGTTATGAACATTTATN",
        "4=1X19=6N10=", "GT", "4S1G29S", "6");

    final HashMap<Integer, String> map = CigarHelper.unrollLegacyCgCigar(r);
    assertTrue(map.size() == 2);
    assertEquals("GACTG", map.get(0));
    assertEquals("TTGAGAGGGANAAAGTTATG      AACATTTATN", map.get(4));



    //GNAAAAGTTATGATTTCAAACGTCCTTTTTAGA 1!02/550-561252316+.-.55.,-+#5245 AS:i:1  NM:i:1  MQ:i:255  GS:Z:AAAA GC:Z:3S2G28S  GQ:Z:,+ XA:i:1  IH:i:1

    //GNAAAAGTTATGATTTCAAACGT      CCTTTTTAGA
    //GNAAA
    //   AAAGTTATGATTTCAAACGT      CCTTTTTAGA
    r = SamTestHelper.getLegacyCGSamRecord("GNAAAAGTTATGATTTCAAACGTCCTTTTTAGA",
        "21=1X1=6N10=", "AAAA", "3S2G28S", ",+");

    final HashMap<Integer, String> map2 = CigarHelper.unrollLegacyCgCigar(r);
    assertTrue(map2.size() == 2);
    assertEquals("GNAAA", map2.get(0));
    assertEquals("AAAGTTATGATTTCAAACGT      CCTTTTTAGA", map2.get(3));

    //39091 131 paolo-bac 370 255 10=6N24=  = 56  -314

    //GAAAGTGGTTACTGTGTGTTGGTGNCTCNGTGGT  /421211303.475100/83,2&.!/./!*5543
    //0123456789012345678901234567890123456789
    //GAAAGTGGTT      ACTGTGTGTTGGTGNCTCNG
    //                                   GTGGT
    //
    r = SamTestHelper.getLegacyCGSamRecord("GAAAGTGGTTACTGTGTGTTGGTGNCTCNGTGGT", "10=6N24=", "GG", "29S1G4S", ",+");

    final HashMap<Integer, String> map3 = CigarHelper.unrollLegacyCgCigar(r);
    assertTrue(map3.size() == 2);
    assertEquals("GAAAGTGGTT      ACTGTGTGTTGGTGNCTCNG", map3.get(0));
    assertEquals("GTGGT", map3.get(35));


    //23663 115 paolo-bac 658 255 10=6N23=  = 308 -349


    //TGTAGTAAAGCTTTCCCCTCTTGAGATTTTANN 45430,-/000...+843/40535.2507.1!!
    //012345678901234567890123456789012345678
    //TGTAGTAAAG      CTTTCCCCTCTTGAGATTTT
    //                                  TTANN
    //
    r = SamTestHelper.getLegacyCGSamRecord("TGTAGTAAAGCTTTCCCCTCTTGAGATTTTANN", "10=6N23=", "AATT", "28S2G3S", "44");

    final HashMap<Integer, String> map4 = CigarHelper.unrollLegacyCgCigar(r);
    assertTrue(map4.size() == 2);
    assertEquals("TGTAGTAAAG      CTTTCCCCTCTTGAGATTTT", map4.get(0));
    assertEquals("AAANN", map4.get(34));
  }

  public void testSuperCigar() {
    final int[] inserts = new int[39];
    CigarHelper.locateInserts("5=2B20=7N1T1R5=2I1=", 1, inserts);
  }

  public void testSoftClipStart() {
    final String cigar = "5S21M1I10M1I6M1D14M1D98M1D79M";
    final int[] inserts = new int[328];
    final int offset = 96;

    CigarHelper.locateInserts(cigar, offset, inserts);

    assertEquals(1, inserts[117]);
    assertEquals(1, inserts[127]);
  }

  public void testInsertsBeforeStartOfTemplate() {
    final String cigar = "2I181M";
    final int[] inserts = new int[232];
    final int offset = 0;
    CigarHelper.locateInserts(cigar, offset, inserts);
    assertEquals(2, inserts[0]);
  }

}
