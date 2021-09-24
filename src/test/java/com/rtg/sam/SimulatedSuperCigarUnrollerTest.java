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
