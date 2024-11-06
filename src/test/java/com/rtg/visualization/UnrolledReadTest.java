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

import com.rtg.mode.DnaUtils;
import com.rtg.sam.BadSuperCigarException;

import htsjdk.samtools.SAMRecord;

import junit.framework.TestCase;

/**
 */
public class UnrolledReadTest extends TestCase {

  public void testUnrolledSuperCigarRead() throws BadSuperCigarException {
    //     0123456789012345678901234567890123456789
    //REF  AACGGATTTAGTGCCCCTGCTCANNNNNNAGTT_ACACA
    //READ AACGG                            _
    //        GGATTTAGTGCCCCTGCTTA      AGTTtACACA

    final SAMRecord rec = SamTestHelper.getSuperCigarSAMRecord("5=21=1X1=6N4=1I5=", "AACGGATTTAGTGCCCCTGCTTAAGTTTACACA", "5=2B18=1X1=6N4=1I5=", "R:", "TT");
    final String template = "AACGGATTTAGTGCCCCTGCTCANNNNNNAGTTACACA";
    final byte[] templateBytes = DnaUtils.encodeString(template);

    final UnrolledRead r = new UnrolledRead(rec, templateBytes, 0);

    assertTrue(r.isInsert(33));
    final String[] res = r.lines(template.length() + 1); //add 1 for only 1 insert
    assertEquals(2, res.length);
    assertEquals("AACGG                                  ", res[0]);
    assertEquals("   GGATTTAGTGCCCCTGCTTA      AGTTtACACA", res[1]);

    final SAMRecord rec2 = SamTestHelper.getSuperCigarSAMRecord("5=21=1X1=6N10=", "AACGGATTTAGTGCCCCTGCTTAAGTTTACACA", "5=2B18=1X1=6N10=", "R:", "T");
    final String template2 = "AACGGATTTAGTGCCCCTGCTCANNNNNNAGTTTACACA";
    final byte[] template2Bytes = DnaUtils.encodeString(template2);

    final UnrolledRead r2 = new UnrolledRead(rec2, template2Bytes, 0);

    assertFalse(r2.isInsert(33));
    final String[] res2 = r2.lines(template.length());
    assertEquals(2, res2.length);
    assertEquals("AACGG                                 ", res2[0]);
    assertEquals("   GGATTTAGTGCCCCTGCTTA      AGTTTACACA", res2[1]);
  }

  public void testLegacyCGCigars() throws BadSuperCigarException {
    //     0123456789012345678901234567890123456789
    // REF GAC_TTGAGAGGGANAAAGTTATGNNNNNNAACATTTATN
    //READ GACtCT
    //          GAGAGGGANAAAGTTATG       AACATTTATN
    //
    //
    final SAMRecord rec = SamTestHelper.getLegacyCGSamRecord("GACTTTGAGAGGGANAAAGTTATGAACATTTATN", "3=1I1X19=6N10=", "GT", "4S1G29S", "6");
    final String template = "GACTTGAGAGGGANAAAGTTATGNNNNNNAACATTTATN";
    final byte[] templateBytes = DnaUtils.encodeString(template);

    final UnrolledRead r = new UnrolledRead(rec, templateBytes, 0);

    assertTrue(r.isInsert(3));

    final String[] res = r.lines(template.length());
    assertEquals(2, res.length);
    assertEquals("GACtG                                  ", res[0]);
    assertEquals("    TTGAGAGGGANAAAGTTATG      AACATTTATN", res[1]);


   // REF GACTTTGAGAGGGANAAAGTTATGAACATTTATN
   //READ GACTT
   //         TGGAGAGGGANAAAGTTATGAACATTTATN
    final SAMRecord rec2 = SamTestHelper.getLegacyCGSamRecord("GACTTTGAGAGGGANAAAGTTATGAACATTTATN", "4=1X19=6N10=", "GT", "4S1G29S", "6");
    final String template2 = "GACTTTGAGAGGGANAAAGTTATGNNNNNNAACATTTATN";
    final byte[] template2Bytes = DnaUtils.encodeString(template2);

    final UnrolledRead r2 = new UnrolledRead(rec2, template2Bytes, 0);

    assertFalse(r2.isInsert(3));
    final String[] res2 = r2.lines(template.length());
    assertEquals(2, res2.length);
    assertEquals("GACTG                                  ", res2[0]);
    assertEquals("    TTGAGAGGGANAAAGTTATG      AACATTTATN", res2[1]);
  }

  public void testCigars() throws BadSuperCigarException {
    //  012345678901234567890123456789012345678901234567890123456789
    // REF AAACCTCTACGCACCCGTCCTCTTTCCGCT_TCTTTCTAGCGGATTTCCCCTAGCTCATC
    //READ AAACCTCTACGCACCCGTCCTCTTTCCGCTGTCTTTCTAGCGGATTTCCCCTAGCTCATC

    final SAMRecord rec = SamTestHelper.getSamRecord("AAACCTCTACGCACCCGTCCTCTTTCCGCTGTCTTTCTAGCGGATTTCCCCTAGCTCATC", "30=1I29=");
    final String template = "AAACCTCTACGCACCCGTCCTCTTTCCGCTTCTTTCTAGCGGATTTCCCCTAGCTCATC";
    final byte[] templateBytes = DnaUtils.encodeString(template);

    final UnrolledRead r = new UnrolledRead(rec, templateBytes, 0);

    assertTrue(r.isInsert(30));
    final String[] res = r.lines(template.length() + 1); //add 1 for insert
    assertEquals(1, res.length);
    assertEquals("AAACCTCTACGCACCCGTCCTCTTTCCGCTgTCTTTCTAGCGGATTTCCCCTAGCTCATC", res[0]);

    final SAMRecord rec2 = SamTestHelper.getSamRecord("AAACCTCTACGCACCCGTCCTCTTTCCGCTTCTTTCTAGCGGATTTCCCCTAGCTCATC", "59=");
    final String template2 = "AAACCTCTACGCACCCGTCCTCTTTCCGCTTCTTTCTAGCGGATTTCCCCTAGCTCATC";
     final UnrolledRead r2 = new UnrolledRead(rec2, DnaUtils.encodeString(template2), 0);

     assertFalse(r2.isInsert(30));
     final String[] res2 = r2.lines(template.length());
     assertEquals(1, res2.length);
     assertEquals(template2, res2[0]);
  }
}
