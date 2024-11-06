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
package com.rtg.variant.util;

import com.rtg.mode.DnaUtils;
import com.rtg.sam.BadSuperCigarException;
import com.rtg.sam.SamUtils;
import com.rtg.variant.VariantAlignmentRecord;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
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
    assertEquals("GAGGCCGAGGCAGGCGGATCGTCAGGAGTT", DnaUtils.bytesToSequenceIncCG(validator.getByteArray()));

    samrec.setCigarString("2=1X7=5N20=");
    samrec.setReadString("GACGCCGAGGCAGGCGGATCGTCAGGAGTT");
    samrec.setAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES, 1);
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "2=1X7=5N20=");
    samrec.setAttribute(SamUtils.CG_READ_DELTA, "C");
    validator.setAlignmentRecord(new VariantAlignmentRecord(samrec));
    validator.parse();
    assertEquals("GACGCCGAGGCAGGCGGATCGTCAGGAGTT", DnaUtils.bytesToSequenceIncCG(validator.getByteArray()));

    samrec.setAlignmentStart(5);
    samrec.setCigarString("4=1D6=5N20=");
    samrec.setReadString("CTGTCATCTTACCTGGGGCCCTCNCTGAGT");
    samrec.setAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES, 1);
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "4=1D6=5N20=");
    samrec.setAttribute(SamUtils.CG_READ_DELTA, null);
    validator.setAlignmentRecord(new VariantAlignmentRecord(samrec));
    validator.setTemplate(DnaUtils.encodeString("TGTTCTGTGCATCTTCCCTTACCTGNGGCCCTCACTGAGTGGGTCCTCCATGGGTGACTGGTGA"));
    validator.parse();
    assertEquals("CTGTCATCTTACCTGNGGCCCTCACTGAGT", DnaUtils.bytesToSequenceIncCG(validator.getByteArray()));

    samrec.setAlignmentStart(5);
    samrec.setCigarString("4=1I5=7N20=");
    samrec.setReadString("CTGTAGCATCACCTGGGGCCCTCNCTGAGT");
    samrec.setAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES, 1);
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "4=1I5=7N20=");
    samrec.setAttribute(SamUtils.CG_READ_DELTA, "A");
    validator.setAlignmentRecord(new VariantAlignmentRecord(samrec));
    validator.parse();
    assertEquals("CTGTAGCATCACCTGNGGCCCTCACTGAGT", DnaUtils.bytesToSequenceIncCG(validator.getByteArray()));

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
    assertEquals("CTTCACAGCGATGGAGAAACTCGGGTGTCTACGTA", DnaUtils.bytesToSequenceIncCG(validator.getByteArray()));

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
    assertEquals("TTTGTGTAGGTCGGATAAGGCGTTCGGATCCGACACG", DnaUtils.bytesToSequenceIncCG(validator.getByteArray()));

    samrec.setCigarString("3=1X21=5N10=");
    samrec.setReadString("tttataggtcggataaggcgttcggatccgacacg");
    samrec.setBaseQualityString("431%68883-56+141663,2.3----45/.,2553");
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "5=2B1X21=5N10=");
    samrec.setAttribute(SamUtils.CG_READ_DELTA, "A");

    validator.setAlignmentRecord(new VariantAlignmentRecord(samrec));
    validator.parse();
    assertEquals("TTTGTATAGGTCGGATAAGGCGTTCGGATCCGACACG", DnaUtils.bytesToSequenceIncCG(validator.getByteArray()));
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
    assertEquals("GATGCCGAGGCAGGCGGATCGTCAGGAGTT", DnaUtils.bytesToSequenceIncCG(validator.getByteArray()));
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
    assertEquals("TTTTTCTTTTTTTTTATTTTTTTTCGTTTTTTTTT", DnaUtils.bytesToSequenceIncCG(validator.getByteArray()));
  }

  public void testOffTemplateSoftClipEnd() throws BadSuperCigarException {
    final SuperCigarUnroller validator = new SuperCigarUnroller();
    final SAMRecord samrec = new SAMRecord(new SAMFileHeader());

    samrec.setAlignmentStart(1);
    samrec.setCigarString("10=7N19=4S");
    samrec.setReadString("AGCCCACACG       TAAATAAGACATCACGATG ATCA".replaceAll(" ", ""));
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "10=7N19=1S2B1=4S");
    samrec.setBaseQualityString("8:::::79:775986<=<<96576767679808");
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR_OVERLAP_QUALITY, "88");
    samrec.setAttribute(SamUtils.CG_READ_DELTA, "AATCA");
//    samrec.setFlags(115);
    validator.setTemplate(DnaUtils.encodeString("AGCCCACACGTTCCCCTTAAATAAGACATCACGATG"));
    validator.setAlignmentRecord(new VariantAlignmentRecord(samrec));
    validator.parse();
    assertEquals("AGCCCACACG TAAATAAGACATCACGATG A GA TCA".replaceAll(" ", ""), DnaUtils.bytesToSequenceIncCG(validator.getByteArray()));
  }

  public void testOffTemplateSoftClipFront() throws BadSuperCigarException {
    final SuperCigarUnroller validator = new SuperCigarUnroller();
    final SAMRecord samrec = new SAMRecord(new SAMFileHeader());

    samrec.setAlignmentStart(1);
    samrec.setCigarString("4S19=7N10=");
    samrec.setReadString("ACTA GTAGCACTACAGAATAAAT       GCACACCCGA".replaceAll(" ", ""));
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "4S1=2B1S19=7N10=");
    samrec.setBaseQualityString("8:::::79:775986<=<<96576767679808");
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR_OVERLAP_QUALITY, "88");
    samrec.setAttribute(SamUtils.CG_READ_DELTA, "ACTAA");
//    samrec.setFlags(115);
    validator.setTemplate(DnaUtils.encodeString("GTAGCACTACAGAATAAATTCCCCTTGCACACCCGA"));
    validator.setAlignmentRecord(new VariantAlignmentRecord(samrec));
    validator.parse();
    assertEquals("ACT AG A GTAGCACTACAGAATAAAT GCACACCCGA".replaceAll(" ", ""), DnaUtils.bytesToSequenceIncCG(validator.getByteArray()));
  }

  public void testReadSoftClipStart() throws Exception {
    final SuperCigarUnroller validator = new SuperCigarUnroller();
    final SAMRecord samrec = new SAMRecord(new SAMFileHeader());

    samrec.setAlignmentStart(5);
    samrec.setCigarString("4S6=7N23=");
    samrec.setReadString("NNNNCACACG       TAAATAAGACATCACGAT GA TCA".replaceAll(" ", ""));
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "4S6=7N20=2B5=");
    samrec.setBaseQualityString("8:::::79:775986<=<<96576767679808");
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR_OVERLAP_QUALITY, "88");
    samrec.setAttribute(SamUtils.CG_READ_DELTA, "NNNN");
//    samrec.setFlags(115);
    validator.setTemplate(DnaUtils.encodeString("TTTTCACACGTTCCCCTTAAATAAGACATCACGATGATCA"));
    validator.setAlignmentRecord(new VariantAlignmentRecord(samrec));
    validator.parse();
    assertEquals("NNNNCACACG TAAATAAGACATCACGAT GA GATCA".replaceAll(" ", ""), DnaUtils.bytesToSequenceIncCG(validator.getByteArray()));
  }
}
