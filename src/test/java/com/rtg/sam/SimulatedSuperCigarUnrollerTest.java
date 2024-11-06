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
package com.rtg.sam;

import com.rtg.mode.DnaUtils;

import junit.framework.TestCase;


/**
 */
public class SimulatedSuperCigarUnrollerTest extends TestCase {


  public void testSimpleMatches() throws Exception {
    final SimulatedSuperCigarUnroller validator = new SimulatedSuperCigarUnroller();


    validator.setTemplateStart(0);
    validator.setTemplate(DnaUtils.encodeString("GAGGCCGAGGGGGGGCAGGCGGATCGTCAGGAGTT"));
    validator.setSdfRead("GAGGCCGAGGCAGGCGGATCGTCAGGAGTT");
    validator.setCigar("10=5N20=", "");
    validator.parse();
    assertEquals("GAGGCCGAGGCAGGCGGATCGTCAGGAGTT", validator.getString());

    validator.setTemplateStart(0);
    validator.setSdfRead("GACGCCGAGGCAGGCGGATCGTCAGGAGTT");
    validator.setCigar("2=1X7=5N20=", "");
    validator.parse();
    assertEquals("GACGCCGAGGCAGGCGGATCGTCAGGAGTT", validator.getString());

    validator.setTemplateStart(4);
    validator.setTemplate(DnaUtils.encodeString("TGTTCTGTGCATCTTCCCTTACCTGNGGCCCTCACTGAGTGGGTCCTCCATGGGTGACTGGTGA"));
    validator.setSdfRead("CTGTCATCTTACCTGNGGCCCTCACTGAGT");
    validator.setCigar("4=1D6=5N20=", "");
    validator.parse();
    assertEquals("CTGTCATCTTACCTGNGGCCCTCACTGAGT", validator.getString());


    /*
    samrec.setAlignmentStart(5);
    samrec.setCigarString("4=1I5=7N20=");
    samrec.setCgReadString("CTGTAGCATCACCTGGGGCCCTCNCTGAGT");
    samrec.setAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES, 1);
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "4=1I5=7N20=");
    samrec.setAttribute(SamUtils.CG_READ_DELTA, "A");
    validator.setSamRecord(samrec);
    validator.parse();
    assertEquals("CTGTAGCATC       ACCTGNGGCCCTCACTGAGT", validator.getString());

     *
     */
    //reverse complement
    /*
    samrec.setCgReadString("CTTCAGCGATGGAGAAACTCGGGTGTCTACGTA");
    samrec.setAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES, 0);
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "5=2B20=6N10=");
    samrec.setAttribute(SamUtils.CG_READ_DELTA, null);
    samrec.setFlags(179);
    validator.setSamRecord(samrec);
    validator.setTemplate(DnaUtils.encodeString("GCTTCAGCGATGGAGAAACTCGGGAAGTCGTGTCTACGTAGAACGTAGTT"));
    validator.parse();
    assertEquals("CTTCACAGCGATGGAGAAACTCGGG      TGTCTACGTA", validator.getString());
     *
     */

    validator.setTemplateStart(1);
    validator.setTemplate(DnaUtils.encodeString("GCTTCAGCGATGGAGAAACTCGGGAAGTCGTGTCTACGTAGAACGTAGTT"));
    validator.setSdfRead("TACGTAGACACCCGAGTTTCTGCATCGCTGTGAAG");
    validator.setReverse(true);
    validator.setCigar("5=2B10=1X9=6N10=", "");
    validator.parse();
    assertEquals("CTTCACAGCGATGGAGAAACTCGGGTGTCTACGTA", validator.getString());


  }

  public void testCgOverlap() throws Exception {
    final SimulatedSuperCigarUnroller validator = new SimulatedSuperCigarUnroller();

    //qual "4316%%68883-56+141663,2.3----45/.,2553"
    validator.setTemplateStart(0);
    validator.setTemplate(DnaUtils.encodeString("tttgtaggtcggataaggcgttcgggggggatccgacacg"));
    validator.setSdfRead("TTTGTGTAGGTCGGATAAGGCGTTCGGATCCGACACG");
    validator.setCigar("5=2B22=5N10=", "");
    validator.parse();
    assertEquals("TTTGTGTAGGTCGGATAAGGCGTTCGGATCCGACACG", validator.getString());

    validator.setTemplateStart(0);
    validator.setTemplate(DnaUtils.encodeString("tttgtaggtcggataaggcgttcgggggggatccgacacg"));
    validator.setSdfRead("TTTGTATAGGTCGGATAAGGCGTTCGGATCCGACACG");
    validator.setCigar("5=2B1X21=5N10=", "");
    validator.parse();
    assertEquals("TTTGTATAGGTCGGATAAGGCGTTCGGATCCGACACG", validator.getString());
  }

  public void testHeapsOfMismatches() throws Exception {
    final SimulatedSuperCigarUnroller validator = new SimulatedSuperCigarUnroller();

    validator.setTemplateStart(0);
    validator.setTemplate(DnaUtils.encodeString("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
    validator.setSdfRead("TTTTTCTTTTTTTTTATTTTTTTTCGTTTTTTTTT");
    validator.setCigar("5X2B10X1=9X5N10X", "");
    validator.parse();
    assertEquals("TTTTTCTTTTTTTTTATTTTTTTTCGTTTTTTTTT", validator.getString());
  }

  public void testOldCigar() {
    try {
      final SimulatedSuperCigarUnroller validator = new SimulatedSuperCigarUnroller();
      validator.setCigar("35=:2B0S5S", null);
      fail();
    } catch (RuntimeException e) {
      assertEquals("old style simulated cigar found: 35=:2B0S5S", e.getMessage());
    }
  }
}
