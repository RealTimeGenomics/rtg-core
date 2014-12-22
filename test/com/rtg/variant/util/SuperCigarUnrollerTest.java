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
package com.rtg.variant.util;

import com.rtg.mode.DnaUtils;
import com.rtg.sam.BadSuperCigarException;
import com.rtg.sam.SamUtils;
import com.rtg.variant.VariantAlignmentRecord;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;

import junit.framework.TestCase;


/**
 */
public class SuperCigarUnrollerTest extends TestCase {


  public void testSimpleMatches() throws Exception {
    final SuperCigarUnroller validator = new SuperCigarUnroller();

    final SAMRecord samrec = new SAMRecord(new SAMFileHeader());

    samrec.setAlignmentStart(1);
    samrec.setCigarString("10M5N20M");
    samrec.setReadString("GAGGCCGAGGCAGGCGGATCGTCAGGAGTT");
    samrec.setAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES, 0);
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "10=5N20=");
    validator.setAlignmentRecord(new VariantAlignmentRecord(samrec));
    validator.setTemplate(DnaUtils.encodeString("GAGGCCGAGGGGGGGCAGGCGGATCGTCAGGAGTT"));
    validator.parse();
    assertEquals("GAGGCCGAGGCAGGCGGATCGTCAGGAGTT", validator.getString());

    samrec.setCigarString("2=1X7=5N20=");
    samrec.setReadString("GACGCCGAGGCAGGCGGATCGTCAGGAGTT");
    samrec.setAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES, 1);
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "2=1X7=5N20=");
    samrec.setAttribute(SamUtils.CG_READ_DELTA, "C");
    validator.setAlignmentRecord(new VariantAlignmentRecord(samrec));
    validator.parse();
    assertEquals("GACGCCGAGGCAGGCGGATCGTCAGGAGTT", validator.getString());

    samrec.setAlignmentStart(5);
    samrec.setCigarString("4=1D6=5N20=");
    samrec.setReadString("CTGTCATCTTACCTGGGGCCCTCNCTGAGT");
    samrec.setAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES, 1);
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "4=1D6=5N20=");
    samrec.setAttribute(SamUtils.CG_READ_DELTA, null);
    validator.setAlignmentRecord(new VariantAlignmentRecord(samrec));
    validator.setTemplate(DnaUtils.encodeString("TGTTCTGTGCATCTTCCCTTACCTGNGGCCCTCACTGAGTGGGTCCTCCATGGGTGACTGGTGA"));
    validator.parse();
    assertEquals("CTGTCATCTTACCTGNGGCCCTCACTGAGT", validator.getString());

    samrec.setAlignmentStart(5);
    samrec.setCigarString("4=1I5=7N20=");
    samrec.setReadString("CTGTAGCATCACCTGGGGCCCTCNCTGAGT");
    samrec.setAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES, 1);
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "4=1I5=7N20=");
    samrec.setAttribute(SamUtils.CG_READ_DELTA, "A");
    validator.setAlignmentRecord(new VariantAlignmentRecord(samrec));
    validator.parse();
    assertEquals("CTGTAGCATCACCTGNGGCCCTCACTGAGT", validator.getString());

    //reverse complement
    samrec.setAlignmentStart(2);
    samrec.setCigarString("23=6N10=");
    samrec.setReadString("CTTCAGCGATGGAGAAACTCGGGTGTCTACGTA");
    samrec.setAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES, 0);
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "5=2B20=6N10=");
    samrec.setAttribute(SamUtils.CG_READ_DELTA, null);
    samrec.setFlags(179);
    validator.setAlignmentRecord(new VariantAlignmentRecord(samrec));
    validator.setTemplate(DnaUtils.encodeString("GCTTCAGCGATGGAGAAACTCGGGAAGTCGTGTCTACGTAGAACGTAGTT"));
    validator.parse();
    assertEquals("CTTCACAGCGATGGAGAAACTCGGGTGTCTACGTA", validator.getString());

  }

  public void testCgOverlap() throws Exception {
    final SuperCigarUnroller validator = new SuperCigarUnroller();
    final SAMRecord samrec = new SAMRecord(new SAMFileHeader());

    samrec.setAlignmentStart(1);
    samrec.setCigarString("25=5N10=");
    samrec.setReadString("tttgtaggtcggataaggcgttcggatccgacacg");
    samrec.setBaseQualityString("431%68883-56+141663,2.3----45/.,2553");
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "5=2B22=5N10=");
    samrec.setAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES, 0);
    samrec.setAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_QUALITY, "6%");
    samrec.setReadNegativeStrandFlag(false);
    samrec.setFlags(67);

    validator.setAlignmentRecord(new VariantAlignmentRecord(samrec));             //   tttgt  aggtcggataaggcgttcgg     atccgacacg
    validator.setTemplate(DnaUtils.encodeString("tttgtaggtcggataaggcgttcgggggggatccgacacg"));
    //qual "4316%%68883-56+141663,2.3----45/.,2553"
    validator.parse();
    assertEquals("TTTGTGTAGGTCGGATAAGGCGTTCGGATCCGACACG", validator.getString());

    samrec.setCigarString("3=1X21=5N10=");
    samrec.setReadString("tttataggtcggataaggcgttcggatccgacacg");
    samrec.setBaseQualityString("431%68883-56+141663,2.3----45/.,2553");
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "5=2B1X21=5N10=");
    samrec.setAttribute(SamUtils.CG_READ_DELTA, "A");

    validator.setAlignmentRecord(new VariantAlignmentRecord(samrec));
    validator.parse();
    assertEquals("TTTGTATAGGTCGGATAAGGCGTTCGGATCCGACACG", validator.getString());
  }

  public void testMismatchFailures() throws Exception {
    final SuperCigarUnroller validator = new SuperCigarUnroller();
    final SAMRecord samrec = new SAMRecord(new SAMFileHeader());

    samrec.setAlignmentStart(1);
    samrec.setCigarString("2=1X7=5N20=");
    samrec.setReadString("GACGCCGAGGCAGGCGGATCGTCAGGAGTT");
    samrec.setAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES, 1);
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "2=1X7=5N20=");
    validator.setAlignmentRecord(new VariantAlignmentRecord(samrec));
    validator.setTemplate(DnaUtils.encodeString("GAGGCCGAGGGGGGGCAGGCGGATCGTCAGGAGTT"));


    samrec.setAttribute(SamUtils.CG_READ_DELTA, "T");
    validator.setAlignmentRecord(new VariantAlignmentRecord(samrec));
    validator.parse();
    assertEquals("GATGCCGAGGCAGGCGGATCGTCAGGAGTT", validator.getString());
  }

  public void testAllMismatches() throws Exception {
    final SuperCigarUnroller validator = new SuperCigarUnroller();
    final SAMRecord samrec = new SAMRecord(new SAMFileHeader());

    samrec.setAlignmentStart(1);
    samrec.setCigarString("2=1X7=5N20=");
    samrec.setReadString("GACGCCGAGGCAGGCGGATCGTCAGGAGTT");
    samrec.setAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES, 1);
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "5X2B10X1=9X5N10X");
    validator.setAlignmentRecord(new VariantAlignmentRecord(samrec));
    validator.setTemplate(DnaUtils.encodeString("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));


    samrec.setAttribute(SamUtils.CG_READ_DELTA, "TTTTTCTTTTTTTTTTTTTTTTTCGTTTTTTTTT");
    validator.setAlignmentRecord(new VariantAlignmentRecord(samrec));
    validator.parse();
    assertEquals("TTTTTCTTTTTTTTTATTTTTTTTCGTTTTTTTTT", validator.getString());
  }

  public void testOffTemplateSoftClipEnd() throws BadSuperCigarException {
    final SuperCigarUnroller validator = new SuperCigarUnroller();
    final SAMRecord samrec = new SAMRecord(new SAMFileHeader());

    samrec.setAlignmentStart(1);
    samrec.setCigarString("10=7N19=4S");
    samrec.setReadString("AGCCCACACG       TAAATAAGACATCACGATG ATCA".replaceAll(" ", ""));
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "10=7N19=1S2B1=4S");
    samrec.setBaseQualityString("8:::::79:775986<=<<96576767679808");
    samrec.setAttribute(SamUtils.CG_OVERLAP_QUALITY, "88");
    samrec.setAttribute(SamUtils.CG_READ_DELTA, "AATCA");
//    samrec.setFlags(115);
    validator.setTemplate(DnaUtils.encodeString("AGCCCACACGTTCCCCTTAAATAAGACATCACGATG"));
    validator.setAlignmentRecord(new VariantAlignmentRecord(samrec));
    validator.parse();
    assertEquals("AGCCCACACG TAAATAAGACATCACGATG A GA TCA".replaceAll(" ", ""), validator.getString());
  }

  public void testOffTemplateSoftClipFront() throws BadSuperCigarException {
    final SuperCigarUnroller validator = new SuperCigarUnroller();
    final SAMRecord samrec = new SAMRecord(new SAMFileHeader());

    samrec.setAlignmentStart(1);
    samrec.setCigarString("4S19=7N10=");
    samrec.setReadString("ACTA GTAGCACTACAGAATAAAT       GCACACCCGA".replaceAll(" ", ""));
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "4S1=2B1S19=7N10=");
    samrec.setBaseQualityString("8:::::79:775986<=<<96576767679808");
    samrec.setAttribute(SamUtils.CG_OVERLAP_QUALITY, "88");
    samrec.setAttribute(SamUtils.CG_READ_DELTA, "ACTAA");
//    samrec.setFlags(115);
    validator.setTemplate(DnaUtils.encodeString("GTAGCACTACAGAATAAATTCCCCTTGCACACCCGA"));
    validator.setAlignmentRecord(new VariantAlignmentRecord(samrec));
    validator.parse();
    assertEquals("ACT AG A GTAGCACTACAGAATAAAT GCACACCCGA".replaceAll(" ", ""), validator.getString());
  }

  public void testReadSoftClipStart() throws Exception {
    final SuperCigarUnroller validator = new SuperCigarUnroller();
    final SAMRecord samrec = new SAMRecord(new SAMFileHeader());

    samrec.setAlignmentStart(5);
    samrec.setCigarString("4S6=7N23=");
    samrec.setReadString("NNNNCACACG       TAAATAAGACATCACGAT GA TCA".replaceAll(" ", ""));
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "4S6=7N20=2B5=");
    samrec.setBaseQualityString("8:::::79:775986<=<<96576767679808");
    samrec.setAttribute(SamUtils.CG_OVERLAP_QUALITY, "88");
    samrec.setAttribute(SamUtils.CG_READ_DELTA, "NNNN");
//    samrec.setFlags(115);
    validator.setTemplate(DnaUtils.encodeString("TTTTCACACGTTCCCCTTAAATAAGACATCACGATGATCA"));
    validator.setAlignmentRecord(new VariantAlignmentRecord(samrec));
    validator.parse();
    assertEquals("NNNNCACACG TAAATAAGACATCACGAT GA GATCA".replaceAll(" ", ""), validator.getString());
  }
}
