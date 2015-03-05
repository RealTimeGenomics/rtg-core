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
    SAMRecord rec = SamTestHelper.getLegacyCGSamRecord("GACTTTGAGAGGGANAAAGTTATGAACATTTATN", "3=1I1X19=6N10=", "GT", "4S1G29S", "6");
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
    SAMRecord rec2 = SamTestHelper.getLegacyCGSamRecord("GACTTTGAGAGGGANAAAGTTATGAACATTTATN", "4=1X19=6N10=", "GT", "4S1G29S", "6");
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
