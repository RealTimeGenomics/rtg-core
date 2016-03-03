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
import com.rtg.reader.FastaUtils;

import htsjdk.samtools.SAMRecord;

import junit.framework.TestCase;


/**
 */
public class SuperCigarValidatorTest extends TestCase {


  public void testSimpleMatches() throws Exception {
    final SuperCigarValidator validator = new SuperCigarValidator(0);

    final SAMRecord samrec = new SAMRecord(null);

    samrec.setAlignmentStart(1);
    samrec.setCigarString("10M5N25M");
    samrec.setReadString("GAGGCCGAGGCAGGCGGATCGTCAGGAGTTAAAAA");
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "10=5N25=");
    samrec.setAttribute(SamUtils.ATTRIBUTE_ALIGNMENT_SCORE, 0);
    samrec.setBaseQualityString("4316%%68883-56+141663,2.3----45/.,2");
    samrec.setFlags(73);
    validator.setData(samrec, DnaUtils.encodeString("GAGGCCGAGG     CAGGCGGATCGTCAGGAGTTAAAAA".replaceAll(" ", "")), FastaUtils.asciiToRawQuality("4316%%68883-56+141663,2.3----45/.,2"));
    validator.setTemplate(DnaUtils.encodeString("GAGGCCGAGGGGGGGCAGGCGGATCGTCAGGAGTTAAAAA"));
    validator.parse();
    assertTrue(validator.getInvalidReason(), validator.isValid());

    samrec.setCigarString("2=1X7=5N25=");
    samrec.setReadString("GACGCCGAGGCAGGCGGATCGTCAGGAGTTAAAAA");
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "2=1X7=5N25=");
    samrec.setAttribute(SamUtils.CG_READ_DELTA, "C");
    validator.setData(samrec, DnaUtils.encodeString("GACGCCGAGG     CAGGCGGATCGTCAGGAGTTAAAAA".replaceAll(" ", "")), FastaUtils.asciiToRawQuality("4316%%68883-56+141663,2.3----45/.,2"));
    validator.parse();
    assertFalse(validator.getInvalidReason(), validator.isValid());

    samrec.setAttribute(SamUtils.ATTRIBUTE_ALIGNMENT_SCORE, null);
    validator.setData(samrec, DnaUtils.encodeString("GACGCCGAGG     CAGGCGGATCGTCAGGAGTTAAAAA".replaceAll(" ", "")), FastaUtils.asciiToRawQuality("4316%%68883-56+141663,2.3----45/.,2"));
    validator.parse();
    assertTrue(validator.getInvalidReason(), validator.isValid());

    samrec.setAlignmentStart(5);
    samrec.setCigarString("4=1D6=5N25=");
    samrec.setReadString("CTGTCATCTTACCTGGGGCCCTCNCTGAGTGGGTC");
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "4=1D6=5N25=");
    samrec.setAttribute(SamUtils.CG_READ_DELTA, null);
    samrec.setAttribute(SamUtils.ATTRIBUTE_ALIGNMENT_SCORE, 2);
    validator.setData(samrec, DnaUtils.encodeString("CTGTCATCTT     ACCTGGGGCCCTCNCTGAGTGGGTC".replaceAll(" ", "")), FastaUtils.asciiToRawQuality("4316%%68883-56+141663,2.3----45/.,2"));
    validator.setTemplate(DnaUtils.encodeString("TGTTCTGTGCATCTTCCCTTACCTGNGGCCCTCACTGAGTGGGTCCTCCATGGGTGACTGGTGA"));
    validator.parse();
    assertTrue(validator.getInvalidReason(), validator.isValid());

    samrec.setAlignmentStart(5);
    samrec.setCigarString("4=1I5=7N25=");
    samrec.setReadString("CTGTAGCATCACCTGGGGCCCTCNCTGAGTGGGTC");
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "4=1I5=7N25=");
    samrec.setAttribute(SamUtils.CG_READ_DELTA, "A");
    samrec.setAttribute(SamUtils.ATTRIBUTE_ALIGNMENT_SCORE, null);
    validator.setData(samrec, DnaUtils.encodeString("CTGTAGCATC     ACCTGGGGCCCTCNCTGAGTGGGTC".replaceAll(" ", "")), FastaUtils.asciiToRawQuality("4316%%68883-56+141663,2.3----45/.,2"));
    validator.parse();
    assertTrue(validator.getInvalidReason(), validator.isValid());

    samrec.setAlignmentStart(5);
    samrec.setCigarString("4=2D6=5N25=");
    samrec.setReadString("CTGTCATCTTACCTGGGGCCCTCNCTGAGTGGGTC");
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "4=2D6=5N25=");
    samrec.setAttribute(SamUtils.CG_READ_DELTA, null);
    samrec.setAttribute(SamUtils.ATTRIBUTE_ALIGNMENT_SCORE, 3);
    validator.setData(samrec, DnaUtils.encodeString("CTGTCATCTT     ACCTGGGGCCCTCNCTGAGTGGGTC".replaceAll(" ", "")), FastaUtils.asciiToRawQuality("4316%%68883-56+141663,2.3----45/.,2"));
    validator.setTemplate(DnaUtils.encodeString("TGTTCTGTGGCATCTTCCCTTACCTGNGGCCCTCACTGAGTGGGTCCTCCATGGGTGACTGGTGA"));
    validator.parse();
    assertTrue(validator.getInvalidReason(), validator.isValid());

    //reverse complement
    samrec.setAlignmentStart(2);
    samrec.setCigarString("23=6N10=");
    samrec.setReadString("CTTCAGCGATGGAGAAACTCGGGTGTCTACGTA");
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "5=2B20=6N10=");
    samrec.setAttribute(SamUtils.CG_READ_DELTA, null);
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR_OVERLAP_QUALITY, "%6");
    samrec.setAttribute(SamUtils.ATTRIBUTE_ALIGNMENT_SCORE, 0);
    samrec.setBaseQualityString("4316%8883-56+141663,2.3----45/.,2");
    samrec.setFlags(179);
    validator.setData(samrec, DnaUtils.encodeString("TACGTAGACA     CCCGAGTTTCTCCATCGCTGTGAAG".replaceAll(" ", "")), FastaUtils.asciiToRawQuality("2,./54----3.2,366141+65-38886%%6134"));
    validator.setTemplate(DnaUtils.encodeString("GCTTCAGCGATGGAGAAACTCGGGAAGTCGTGTCTACGTAGAACGTAGTT"));
    validator.parse();
    assertTrue(validator.getInvalidReason(), validator.isValid());

  }

  public void testCGOverlapWithDeletion() throws Exception {
    //check that it is OK to not provide XQ if the overlap is deleted from the template
    final SuperCigarValidator validator = new SuperCigarValidator(0);

    final SAMRecord samrec = new SAMRecord(null);

    samrec.setAlignmentStart(1);
    samrec.setCigarString("29=");
    samrec.setReadString("AGGCAGGTAGATCATGAGGTGAAGAGATC");
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "19=2B2D10=");
    samrec.setAttribute(SamUtils.ATTRIBUTE_ALIGNMENT_SCORE, 3);
    samrec.setBaseQualityString("/////////////////////////////");
    samrec.setFlags(179);
    final byte[] sdfRead = DnaUtils.encodeString("GATCTCTTCACCTCATGATCTACCTGCCT".replaceAll(" ", ""));
    final byte[] sdfQualities = FastaUtils.asciiToRawQuality("/////////////////////////////");
    validator.setData(samrec, sdfRead, sdfQualities);
    validator.setTemplate(DnaUtils.encodeString("AGGCAGGTAGATCATGAGGTGAAGAGATC"));
    validator.parse();
    assertTrue(validator.getInvalidReason(), validator.isValid());

    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "19=2B2N10=");
    samrec.setAttribute(SamUtils.ATTRIBUTE_ALIGNMENT_SCORE, 0);
    validator.setData(samrec, sdfRead, sdfQualities);
    validator.parse();
    assertTrue(validator.getInvalidReason(), validator.isValid());

    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "19=2B2H10=");
    samrec.setAttribute(SamUtils.ATTRIBUTE_ALIGNMENT_SCORE, 0);
    validator.setData(samrec, sdfRead, sdfQualities);
    validator.parse();
    assertTrue(validator.getInvalidReason(), validator.isValid());


    samrec.setCigarString("19=2D8=");
    samrec.setReadString("TGGCAGGTAGATCATGAGGAAGAGATC"); //<- this doesn't seem to be checked by anything
    samrec.setBaseQualityString("/////////////////////////////");
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "19=2B1X1=2D8=");
    samrec.setAttribute(SamUtils.ATTRIBUTE_ALIGNMENT_SCORE, 4);
    samrec.setAttribute(SamUtils.CG_READ_DELTA, "T");
    validator.setData(samrec, sdfRead, sdfQualities);
    validator.parse();
    assertFalse(validator.isValid());
    assertTrue(validator.getInvalidReason(), validator.getInvalidReason().contains("Overlap described but no XQ field present in SAM record"));

    samrec.setReadString("AGGCAGGTAGATCATGAGGAAGAGATC");
    samrec.setBaseQualityString("/////////////////////////////");
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "19=2B1X1=3X3=2X");
    samrec.setAttribute(SamUtils.CG_READ_DELTA, "TAAGTC");
    samrec.setAttribute(SamUtils.ATTRIBUTE_ALIGNMENT_SCORE, 6);
    validator.setData(samrec, sdfRead, sdfQualities);
    validator.parse();
    assertFalse(validator.isValid());
    assertTrue(validator.getInvalidReason(), validator.getInvalidReason().contains("Overlap described but no XQ field present in SAM record"));

  }

    public void testCGOverlapWithDeletion2() throws Exception {
      //check that it is OK to not provide XQ if the overlap is deleted from the template
      final SuperCigarValidator validator = new SuperCigarValidator(0);

      final SAMRecord samrec = new SAMRecord(null);
//    GGGCCTGCAC
//              DDD
//               BB
//               TGGCCAAGGAGCTGTGTGA
//    GGGCCTGCACCTGGCCAAGGAGCTGTGTGA
//
      samrec.setAlignmentStart(1);
      samrec.setCigarString("10=1D19=");
      samrec.setReadString("GGGCCTGCACTGGCCAAGGAGCTGTGTGA");
      samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "10=3D2B19=");
      samrec.setAttribute(SamUtils.ATTRIBUTE_ALIGNMENT_SCORE, 4);
      samrec.setBaseQualityString("/////////////////////////////");
      samrec.setFlags(131);
      final byte[] sdfRead = DnaUtils.encodeString("GGGCCTGCACTGGCCAAGGAGCTGTGTGA".replaceAll(" ", ""));
      final byte[] sdfQualities = FastaUtils.asciiToRawQuality("/////////////////////////////");
      validator.setData(samrec, sdfRead, sdfQualities);
      validator.setTemplate(DnaUtils.encodeString("GGGCCTGCACCTGGCCAAGGAGCTGTGTGA"));
      validator.parse();
      assertTrue(validator.getInvalidReason(), validator.isValid());

      //theoretical alignment that probably isn't handled:
//    ACG
//       DDD
//          G
//         BB
//         C
//          D
//           TACGTACGTACGT
//    ACGTACGTACGTACGTACGT
//    The overlap actually has a match on either side of it, however no template position is repeated
//    in a match or mismatch so would not result in a flattened read needing an XQ field.
    }

  public void testCgOverlap() throws Exception {
    final SuperCigarValidator validator = new SuperCigarValidator(0);
    final SAMRecord samrec = new SAMRecord(null);

    samrec.setAlignmentStart(1);
    samrec.setCigarString("25=5N10=");
    samrec.setReadString("tttgtaggtcggataaggcgttcatccgacacg");
    samrec.setBaseQualityString("4316%8883-56+141663,2.3----45/.,2");
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "5=2B20=5N10=");
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR_OVERLAP_QUALITY, "%6");
    samrec.setReadNegativeStrandFlag(false);
    samrec.setFlags(67);
    validator.setData(samrec, DnaUtils.encodeString("tttgtgtaggtcggataaggcgttc     atccgacacg".replaceAll(" ", "")), FastaUtils.asciiToRawQuality("4316%%68883-56+141663,2.3----45/.,2"));
    validator.setTemplate(DnaUtils.encodeString("tttgtaggtcggataaggcgttcgggggatccgacacg"));
    validator.parse();
    assertTrue(validator.getInvalidReason(), validator.isValid());

    samrec.setCigarString("3=1X21=5N10=");
    samrec.setReadString("tttataggtcggataaggcgttcatccgacacg");
    samrec.setBaseQualityString("4316%8883-56+141663,2.3----45/.,2");
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "5=2B1X19=5N10=");
    samrec.setAttribute(SamUtils.CG_READ_DELTA, "A");

    validator.setData(samrec, DnaUtils.encodeString("tttgtataggtcggataaggcgttc     atccgacacg".replaceAll(" ", "")), FastaUtils.asciiToRawQuality("4316%%68883-56+141663,2.3----45/.,2"));
    validator.parse();
    assertTrue(validator.getInvalidReason(), validator.isValid());
  }

  public void testMismatchFailures() throws Exception {
    final SuperCigarValidator validator = new SuperCigarValidator(0);
    final SAMRecord samrec = new SAMRecord(null);

    samrec.setAlignmentStart(1);
    samrec.setCigarString("2=1X13=5N20=");
    samrec.setReadString("GACGCCGAGGAAAAACAGGCGGATCGTCAGGAGTT");
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "2=1X13=5N20=");
    samrec.setBaseQualityString("4316%%68883-56+141663,2.3----45/.,2");
    samrec.setFlags(73);
    validator.setData(samrec, DnaUtils.encodeString("GACGCCGAGG     AAAAACAGGCGGATCGTCAGGAGTT".replaceAll(" ", "")), FastaUtils.asciiToRawQuality("4316%%68883-56+141663,2.3----45/.,2"));
    validator.setTemplate(DnaUtils.encodeString("GAGGCCGAGGGGGGGAAAAACAGGCGGATCGTCAGGAGTT"));

    try {
      validator.parse();
      assertFalse(validator.isValid());
      assertTrue(validator.getInvalidReason(), validator.getInvalidReason().contains("Read delta (" + SamUtils.CG_READ_DELTA + ") too short, "));
    } catch (final AssertionError e) {
      assertEquals("readDelta.len=0 but should be 1", e.getMessage());
    }

    samrec.setAttribute(SamUtils.CG_READ_DELTA, "T");
    validator.setData(samrec, DnaUtils.encodeString("GACGCCGAGG     AAAAACAGGCGGATCGTCAGGAGTT".replaceAll(" ", "")), FastaUtils.asciiToRawQuality("4316%%68883-56+141663,2.3----45/.,2"));
    validator.parse();
    assertFalse(validator.isValid());
    assertTrue(validator.getInvalidReason(), validator.getInvalidReason().contains(SamUtils.CG_READ_DELTA + " value: T does not match read value: C"));

    validator.setData(samrec, DnaUtils.encodeString("GAGGCCGAGG     AAAAACAGGCGGATCGTCAGGAGTT".replaceAll(" ", "")), FastaUtils.asciiToRawQuality("4316%%68883-56+141663,2.3----45/.,2"));
    validator.parse();
    assertFalse(validator.isValid());
    assertTrue(validator.getInvalidReason(), validator.getInvalidReason().contains("Expected mismatch"));

    samrec.setAlignmentStart(5);
    samrec.setCigarString("4=1I5=7N25=");
    samrec.setReadString("CTGTGGCATCGGGGGACCTGGGGCCCTCNCTGAGT");
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "4=1I5=7N25=");
    samrec.setAttribute(SamUtils.CG_READ_DELTA, "A");
    validator.setData(samrec, DnaUtils.encodeString("CTGTGGCATC     GGGGGACCTGGGGCCCTCNCTGAGT".replaceAll(" ", "")), FastaUtils.asciiToRawQuality("4316%%68883-56+141663,2.3----45/.,2"));
    validator.setTemplate(DnaUtils.encodeString("TGTTCTGTG CATCTTCCCTTGGGGGACCTGNGGCCCTCACTGAGTGGGTCCTCCATGGGTGACTGGTGA"));
    validator.parse();
    assertFalse(validator.getInvalidReason(), validator.isValid());
    assertTrue(validator.getInvalidReason(), validator.getInvalidReason().contains("SDF read insert: G does not match SAM " + SamUtils.CG_READ_DELTA + ": A,"));

    samrec.setAlignmentStart(2);
    samrec.setCigarString("23=6N10=");
    samrec.setReadString("CTTCAGCGATGGAGAAACTCGGGTGTCTACGTA");
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "5=2B20=6N10=");
    samrec.setAttribute(SamUtils.CG_READ_DELTA, null);
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR_OVERLAP_QUALITY, "%6");
    samrec.setBaseQualityString("4316%8883-56+141663,2.3----45/.,2");
    samrec.setFlags(179);
    validator.setData(samrec, DnaUtils.encodeString("TACGTAGACA     CCCGAGTGTCTCCATCGCTGTGAAG".replaceAll(" ", "")), FastaUtils.asciiToRawQuality("2,./54----3.2,366141+65-38886%%6134"));
    validator.setTemplate(DnaUtils.encodeString("GCTTCAGCGATGGAGAAACTCGGGAAGTCGTGTCTACGTAGAACGTAGTT"));
    validator.parse();
    assertFalse(validator.getInvalidReason(), validator.isValid());
    assertTrue(validator.getInvalidReason(), validator.getInvalidReason().contains("Expected match, SDF read=C, template=A,"));

    samrec.setAlignmentStart(1);
    samrec.setCigarString("10M5N25M");
    samrec.setReadString("GAGGCCGAGGCAGGCGGATCGTCAGGAGTTAAAAA");
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "10=5N25=");
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR_OVERLAP_QUALITY, null);
    samrec.setBaseQualityString("4316%668883-56+141663,2.3----45/.,2");
    samrec.setFlags(73);
    validator.setData(samrec, DnaUtils.encodeString("GAGGCCGAGG     CAGGCGGATCGTCAGGAGTTAAAAA".replaceAll(" ", "")), FastaUtils.asciiToRawQuality("4316%%68883-56+141663,2.3----45/.,2"));
    validator.setTemplate(DnaUtils.encodeString("GAGGCCGAGGGGGGGCAGGCGGATCGTCAGGAGTTAAAAA"));
    validator.parse();
    assertFalse(validator.getInvalidReason(), validator.isValid());
    assertTrue(validator.getInvalidReason(), validator.getInvalidReason().contains("SDF and SAM qualities don't match,"));

    samrec.setAlignmentStart(1);
    samrec.setCigarString("25=5N10=");
    samrec.setReadString("tttgtaggtcggataaggcgttcatccgacacg");
    samrec.setBaseQualityString("4316%8883-56+141663,2.3----45/.,2");
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "5=2B20=5N10=");
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR_OVERLAP_QUALITY, "6");
    samrec.setReadNegativeStrandFlag(false);
    samrec.setFlags(67);
    try {
      validator.setData(samrec, DnaUtils.encodeString("tttgtgtaggtcggataaggcgttc     atccgacacg".replaceAll(" ", "")), FastaUtils.asciiToRawQuality("4316%%68883-56+141663,2.3----45/.,2"));
      fail();
    } catch (final BadSuperCigarException iae) {
      assertTrue(iae.getMessage(), iae.getMessage().contains("SAM record qualities plus XQ not expected length. Was: 34 expected: 35"));
    }

    samrec.setAlignmentStart(1);
    samrec.setCigarString("25=5N10=");
    samrec.setReadString("tttgtaggtcggataaggcgttcatccgacacg");
    samrec.setBaseQualityString("4316%%68883-56+141663,2.3----45/.,2");
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "5=2B20=5N10=");
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR_OVERLAP_QUALITY, null);
    samrec.setReadNegativeStrandFlag(false);
    samrec.setFlags(67);
    validator.setTemplate(DnaUtils.encodeString("tttgtaggtcggataaggcgttcgggggatccgacacg"));
    validator.setData(samrec, DnaUtils.encodeString("tttgtgtaggtcggataaggcgttc     atccgacacg".replaceAll(" ", "")), FastaUtils.asciiToRawQuality("4316%%68883-56+141663,2.3----45/.,2"));
    validator.parse();

    assertFalse(validator.getInvalidReason(), validator.isValid());
    assertTrue(validator.getInvalidReason(), validator.getInvalidReason().contains("Overlap described but no " + SamUtils.CG_SUPER_CIGAR_OVERLAP_QUALITY + " field present"));
  }

  public void testQualities() throws Exception {

    //first, non-rc
    /*237726  67      simulatedSequence1      3       255     22=6N10=        =       170     167     TGCCCCCCTGAGAATGAATGTTGGACGAAATA        )*N\S\7@*`[4DRA8VKE-JF:KP0<D:/"K        AS:i:0  NM:i:0  MQ:i:255        XU:Z:5=3B20=6N10=
                                                                                                    TGCCCCCCCCCTGAGAATGAATGTTGGACGAAATA
                                                                                                    )*N\SV55\7@*`[4DRA8VKE-JF:KP0<D:/"K


           XQ:Z:V55        XA:i:1  IH:i:1*/

    final SuperCigarValidator validator = new SuperCigarValidator(0);
    final SAMRecord samrec = new SAMRecord(null);

    samrec.setAlignmentStart(3);
    samrec.setCigarString("22=6N10=");
    samrec.setReadString("tgcccccctgagaatgaatgttggacgaaata");
    samrec.setBaseQualityString(")*N!S!7@*`[4DRA8VKE-JF:KP0<D:/!K");
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "5=3B20=6N10=");
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR_OVERLAP_QUALITY, "V55");
    samrec.setFlags(67);

    validator.setData(samrec, DnaUtils.encodeString("tgccccccccctgagaatgaatgtt     ggacgaaata".replaceAll(" ", "")), FastaUtils.asciiToRawQuality(")*N!SV55!7@*`[4DRA8VKE-JF:KP0<D:/!K"));             //   tttgt  aggtcggataaggcgttcgg     atccgacacg
    validator.setTemplate(DnaUtils.encodeString("attgcccccctgagaatgaatgttatgtacggacgaaatatgtaaccata"));

    validator.parse();
    assertTrue(validator.getInvalidReason(), validator.isValid());

    //second, non-rc
    /*184100  131     simulatedSequence1      191     255     10=6N4=1X18=    =       11      -180    AGCTTCTATAGCGGAATTGAGCGGAACCGCACG       YTD$B1L!`_<"L'V8W=72T#YU]K@,#KUA>       AS:i:1  NM:i:1  MQ:i:255        XU:Z:10=6N4=1X15=2B5=   XQ:Z:%* XR:Z:A  XA:i:1  IH:i:1
                                                                                                    AGCTTCTATAGCGGAATTGAGCGGAACCGCGCACG
                                                                                                    YTD$B1L!`_<"L'V8W=72T#YU]K@,%*#KUA>
     */
    samrec.setAlignmentStart(2);
    samrec.setCigarString("10=6N4=1X18=");
    samrec.setReadString("AGCTTCTATAGCGGAATTGAGCGGAACCGCACG");
    samrec.setBaseQualityString("YTD$B1L!`_<!L'V8W=72T#YU]K@,#KUA>");
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "10=6N4=1X15=2B5=");
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR_OVERLAP_QUALITY, "%*");
    samrec.setAttribute(SamUtils.CG_READ_DELTA, "A");
    samrec.setFlags(131);
    validator.setData(samrec, DnaUtils.encodeString("AGCTTCTATA     GCGGAATTGAGCGGAACCGCGCACG".replaceAll(" ", "")), FastaUtils.asciiToRawQuality("YTD$B1L!`_<!L'V8W=72T#YU]K@,%*#KUA>"));       //    YTD$B1L!`_      <!L'V8W=72T#YU]K@,%*#KUA>
    validator.setTemplate(DnaUtils.encodeString("TAGCTTCTATAGGGGGCGCGGTATTGAGCGGAACCGCACGTGCTATTTTCC"));
    validator.parse();
    assertTrue(validator.getInvalidReason(), validator.isValid());

    //first, rc
    /*137629  115     simulatedSequence1      195     255     10=6N9=1I13=    =       13      -181    TCTATAGGGGTATTGAGCGAGAACCGCACGTGC       ^#R"E\`,),UQANL6J/J"G/P'^;<RIX4O$       AS:i:2  NM:i:1  MQ:i:255        XU:Z:10=6N9=1I10=2B5=   XQ:Z:<O XR:Z:A  XA:i:3  IH:i:1
                                                                                                    GCACGCGTGCGGTTCTCGCTCAATACCCCTATAGA
                                                                                                    $O4XIO<R<;^'P/G"J/J6LNAQU,),`\E"R#^
     */
    samrec.setAlignmentStart(6);
    samrec.setCigarString("10=6N9=1I13=");
    samrec.setReadString("TCTATAGGGGTATTGAGCGAGAACCGCACGTGC");
    samrec.setBaseQualityString("^#R!E!`,),UQANL6J/J!G/P'^;<RIX4O$");
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "10=6N9=1I10=2B5=");
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR_OVERLAP_QUALITY, "<O");
    samrec.setAttribute(SamUtils.CG_READ_DELTA, "A");
    samrec.setFlags(115);
    validator.setData(samrec, DnaUtils.encodeString("GCACGCGTGCGGTTCTCGCTCAATA     CCCCTATAGA".replaceAll(" ", "")), FastaUtils.asciiToRawQuality("$O4XIO<R<;^'P/G!J/J6LNAQU,),`!E!R#^"));       //    YTD$B1L!`_      <!L'V8W=72T#YU]K@,%*#KUA>
    validator.setTemplate(DnaUtils.encodeString("TAGCTTCTATAGGGGGCGCGGTATTGAGCGGAACCGCACGTGCTATTTTCC"));
    validator.parse();
    assertTrue(validator.getInvalidReason(), validator.isValid());



    //second, rc
    /*137629  179     simulatedSequence1      13      255     9=1X13=6N10=    =       195     181     AGAATGAATCTTATGTACGGACGGTAACCATAA       ^#R"E,),UQANL6J/J"G/P'^;<R<OIX4O$       AS:i:1  NM:i:1  MQ:i:255        XU:Z:5=2B6=1X13=6N10=   XQ:Z:\` XR:Z:C  XA:i:3  IH:i:1
                                                                                                    TTATGGTTACCGTCCGTACATAAGATTCATATTCT
                                                                                                    $O4XIO<R<;^'P/G"J/J6LNAQU,),`\E"R#^
                                                                                                    */
    samrec.setAlignmentStart(2);
    samrec.setCigarString("9=1X13=6N10=");
    samrec.setReadString("AGAATGAATCTTATGTACGGACGGTAACCATAA");
    samrec.setBaseQualityString("^#R!E,),UQANL6J/J!G/P'^;<R<OIX4O$");
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "5=2B6=1X13=6N10=");
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR_OVERLAP_QUALITY, "!`");
    samrec.setAttribute(SamUtils.CG_READ_DELTA, "C");
    samrec.setFlags(179);
    validator.setData(samrec, DnaUtils.encodeString("TTATGGTTAC     CGTCCGTACATAAGATTCATATTCT".replaceAll(" ", "")), FastaUtils.asciiToRawQuality("$O4XIO<R<;^'P/G!J/J6LNAQU,),`!E!R#^"));       //    YTD$B1L!`_      <!L'V8W=72T#YU]K@,%*#KUA>
    validator.setTemplate(DnaUtils.encodeString("GAGAATGAATGTTATGTACGGACGAAATATGTAACCATAACACC"));
    validator.parse();
    assertTrue(validator.getInvalidReason(), validator.isValid());
  }

  public void test4Gap() throws Exception {
                                                            //                                               CCATTCAGTTGGAGACGTTGTGGACCTGACGCCTCTGCTCTTGCAAGTCAGGACAT
    //24106       67      paolo-bac       420     255     16=4I1=1I2=4N10=        paolo-bac       735     315     CAGTTGGAGACGTTGTGNATGTGNACGCCTCTGC      213.3/22..103350/!2,2+/!14/-+-4//5      AS:i:6  NM:i:5  MQ:i:255        XU:Z:5=1B12=4I1=1I2=4N10=       XQ:Z:1  XR:Z:GNATT      XA:i:10 IH:i:1
    final SuperCigarValidator validator = new SuperCigarValidator(0);
    final SAMRecord samrec = new SAMRecord(null);

    samrec.setAlignmentStart(6);
    samrec.setCigarString("16=4I1=1I2=4N10=");
    samrec.setReadString("CAGTTGGAGACGTTGTGNATGTGNACGCCTCTGC");
    samrec.setBaseQualityString("213.3/22..103350/!2,2+/!14/-+-4//5");
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "5=1B12=4I1=1I2=4N10=");
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR_OVERLAP_QUALITY, "1");
    samrec.setAttribute(SamUtils.CG_READ_DELTA, "GNATT");
    samrec.setAttribute(SamUtils.ATTRIBUTE_ALIGNMENT_SCORE, 7);
    samrec.setFlags(67);
                                                   // CAGTT GGAGACGTTGTGNATG    T GN   ACGCCTCTGC
    validator.setData(samrec, DnaUtils.encodeString("CAGTTTGGAGACGTTGTGNATGTGN     ACGCCTCTGC".replaceAll(" ", "")), FastaUtils.asciiToRawQuality("213.31/22..103350/!2,2+/!14/-+-4//5"));             //   tttgt  aggtcggataaggcgttcgg     atccgacacg
    validator.setTemplate(DnaUtils.encodeString("CCATTCAGTTGGAGACGTTGTGGACCTGACGCCTCTGCTCTTGCAAGTCAGGACAT"));

    validator.parse();
    assertTrue(validator.getInvalidReason(), validator.isValid());

    samrec.setAttribute(SamUtils.ATTRIBUTE_ALIGNMENT_SCORE, 6);
    validator.setData(samrec, DnaUtils.encodeString("CAGTTTGGAGACGTTGTGNATGTGN     ACGCCTCTGC".replaceAll(" ", "")), FastaUtils.asciiToRawQuality("213.31/22..103350/!2,2+/!14/-+-4//5"));             //   tttgt  aggtcggataaggcgttcgg     atccgacacg
    validator.setTemplate(DnaUtils.encodeString("CCATTCAGTTGGAGACGTTGTGGACCTGACGCCTCTGCTCTTGCAAGTCAGGACAT"));

    validator.parse();
    assertFalse(validator.getInvalidReason(), validator.isValid());
    assertEquals("Super cigar alignment score was 7, but AS attribute was 6", validator.getInvalidReason());
  }

  public void testDegenerate() throws Exception {
    final SuperCigarValidator validator = new SuperCigarValidator(0);
    final SAMRecord samrec = new SAMRecord(null);

    samrec.setAlignmentStart(1);
    samrec.setCigarString("2=1X13=5N20=");
    samrec.setReadString("GACGCCGAGGAAAAACAGGCGGATCGTCAGGAGTT");
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "2=1X13=5N20=");
    samrec.setBaseQualityString("4316%%68883-56+141663,2.3----45/.,2");
    samrec.setFlags(73);
    validator.setData(samrec, DnaUtils.encodeString("                                   "), FastaUtils.asciiToRawQuality("4316%%68883-56+141663,2.3----45/.,2"));
    validator.setTemplate(DnaUtils.encodeString("GAGGCCGAGGGGGGGAAAAACAGGCGGATCGTCAGGAGTT"));

    validator.parse();
    assertFalse(validator.isValid());

    samrec.setFlags(115);
    validator.setData(samrec, DnaUtils.encodeString("                                   "), FastaUtils.asciiToRawQuality("4316%%68883-56+141663,2.3----45/.,2"));

    validator.parse();
    assertFalse(validator.isValid());
  }

  public void testSoftClip() throws Exception {
    final SuperCigarValidator validator = new SuperCigarValidator(0);
    final SAMRecord samrec = new SAMRecord(null);

    samrec.setAlignmentStart(1);
    samrec.setCigarString("3S2=18=5N9=1X");
    samrec.setReadString("AGCCCACACGTAAATAAGACATCACGATGATCA");
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "3S2=2B20=5N9=1X");
    samrec.setBaseQualityString("4316%8883-56+141663,2.3----45/.,2");
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR_OVERLAP_QUALITY, "%6");
    samrec.setAttribute(SamUtils.CG_READ_DELTA, "AGCA");
    samrec.setAttribute(SamUtils.ATTRIBUTE_ALIGNMENT_SCORE, 4);
    samrec.setFlags(67);
    validator.setData(samrec, DnaUtils.encodeString("AGCCCCCACACGTAAATAAGACATC     ACGATGATCA".replaceAll(" ", "")), FastaUtils.asciiToRawQuality("4316%%68883-56+141663,2.3----45/.,2"));
    validator.setTemplate(DnaUtils.encodeString("CCACACGTAAATAAGACATCGGGGGACGATGATCG"));

    validator.parse();
    assertTrue(validator.getInvalidReason(), validator.isValid());


    //Alignment mismatch 3617884      153     chr18   1       255     4S19=7N2=1X7=   *       *       *       AAAACCCTAACCCTAACCCTAACCCCAACCCTA       998140-,7::;26;;.39'(2347-88989+7       AS:i:2  NM:i:1  XU:Z:5=2B20=7N2=1X7=    XR:Z:C  XQ:Z:42 IH:i:1
    samrec.setAlignmentStart(1);
    samrec.setCigarString("4S19=7N2=1X7=");
    samrec.setReadString("AAAACCCTAACCCTAACCCTAACCCCAACCCTA");
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "4S1=2B1S19=7N2=1X7=");
    samrec.setBaseQualityString("998140-,7::;26;;.39'(2347-88989+7");
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR_OVERLAP_QUALITY, "42");
    samrec.setAttribute(SamUtils.CG_READ_DELTA, "AAAAAC");
    samrec.setAttribute(SamUtils.ATTRIBUTE_ALIGNMENT_SCORE, 6);
    samrec.setFlags(67);
    validator.setData(samrec, DnaUtils.encodeString("AAAACACCCTAACCCTAACCCTAAC     CCCAACCCTA".replaceAll(" ", "")), FastaUtils.asciiToRawQuality("99814420-,7::;26;;.39'(2347-88989+7"));
    validator.setTemplate(DnaUtils.encodeString("CCCTAACCCTAACCCTAACCCTTACCCCTAACCCTA"));

    validator.parse();
    assertTrue(validator.getInvalidReason(), validator.isValid());

    //153
    //TAGGGTTGGG     GTTAGGGTTAGGGTTAGGGTGTTTT"), DnaUtils.fastqToPhred("7+98988-7432('93.;;62;::7,-02441899"));


    samrec.setAlignmentStart(1);
    samrec.setCigarString("10=7N19=4S");
    samrec.setReadString("AGCCCACACGTAAATAAGACATCACGATGATCA");
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "10=7N19=1S2B1=4S");
    samrec.setBaseQualityString("8:::::79:775986<=<<96576767679808");
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR_OVERLAP_QUALITY, "88");
    samrec.setAttribute(SamUtils.CG_READ_DELTA, "AATCA");
    samrec.setAttribute(SamUtils.ATTRIBUTE_ALIGNMENT_SCORE, 5);
    samrec.setFlags(139);
    validator.setData(samrec, DnaUtils.encodeString("AGCCCACACG     TAAATAAGACATCACGATGAGATCA".replaceAll(" ", "")), FastaUtils.asciiToRawQuality("8:::::79:775986<=<<9657676768879808"));
    validator.setTemplate(DnaUtils.encodeString("AGCCCACACGTTCCCCTTAAATAAGACATCACGATG"));

    validator.parse();
    assertTrue(validator.getInvalidReason(), validator.isValid());


    samrec.setFlags(115);
    validator.setData(samrec, DnaUtils.encodeString("TGATCTCATCGTGATGTCTTATTTA     CGTGTGGGCT".replaceAll(" ", "")), FastaUtils.asciiToRawQuality("8089788676767569<<=<689577:97:::::8"));

    validator.parse();
    assertTrue(validator.getInvalidReason(), validator.isValid());
  }

  public void testUnknowns() throws Exception {
    final SuperCigarValidator validator = new SuperCigarValidator(0);
    final SAMRecord samrec = new SAMRecord(null);

    samrec.setAlignmentStart(1);
    samrec.setCigarString("8=1X14=5N10=");
    samrec.setReadString("AGCCCCCAAACGTAAATAAGACATCACGATGATCA");
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "7=1R1X1T15=5N10=");
    samrec.setBaseQualityString("4316%%68883-56+141663,2.3----45/.,2");
    samrec.setAttribute(SamUtils.CG_READ_DELTA, "A");
    samrec.setAttribute(SamUtils.ATTRIBUTE_ALIGNMENT_SCORE, 1);
    samrec.setFlags(67);
    validator.setData(samrec, DnaUtils.encodeString("AGCCCCCAAACGTAAATAAGACATC     ACGATGATCA".replaceAll(" ", "")), FastaUtils.asciiToRawQuality("4316%%68883-56+141663,2.3----45/.,2"));
    validator.setTemplate(DnaUtils.encodeString("AGCCCCCACACGTAAATAAGACATCTTTTTACGATGATCA"));

    validator.parse();
    assertFalse(validator.getInvalidReason(), validator.isValid());
    assertEquals("Expected unknown on read, was A", validator.getInvalidReason());

    samrec.setReadString("AGCCCCCNTACGTAAATAAGACATCACGATGATCA");
    validator.setData(samrec, DnaUtils.encodeString("AGCCCCCNAACGTAAATAAGACATC     ACGATGATCA".replaceAll(" ", "")), FastaUtils.asciiToRawQuality("4316%%68883-56+141663,2.3----45/.,2"));
    validator.setTemplate(DnaUtils.encodeString("AGCCCCCACACGTAAATAAGACATCTTTTTACGATGATCA"));

    validator.parse();
    assertFalse(validator.getInvalidReason(), validator.isValid());
    assertTrue(validator.getInvalidReason(), validator.getInvalidReason().startsWith("Read delta (XR) too short"));

    samrec.setAttribute(SamUtils.CG_READ_DELTA, "TA");
    validator.setData(samrec, DnaUtils.encodeString("AGCCCCCNTACGTAAATAAGACATC     ACGATGATCA".replaceAll(" ", "")), FastaUtils.asciiToRawQuality("4316%%68883-56+141663,2.3----45/.,2"));
    validator.setTemplate(DnaUtils.encodeString("AGCCCCCACACGTAAATAAGACATCTTTTTACGATGATCA"));

    validator.parse();
    assertFalse(validator.getInvalidReason(), validator.isValid());
    assertEquals("Expected unknown on template, was A", validator.getInvalidReason());

    validator.setTemplate(DnaUtils.encodeString("AGCCCCCACNCGTAAATAAGACATCTTTTTACGATGATCA"));

    validator.parse();
    assertTrue(validator.getInvalidReason(), validator.isValid());
  }

  public void testUnknownPenalty() throws Exception {
    final SuperCigarValidator validator = new SuperCigarValidator(1);
    final SAMRecord samrec = new SAMRecord(null);

    samrec.setAlignmentStart(1);
    samrec.setCigarString("8=1X14=5N10=");
    samrec.setReadString("AGCCCCCNTACGTAAATAAGACATCACGATGATCA");
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "7=1R1X1T15=5N10=");
    samrec.setBaseQualityString("4316%%68883-56+141663,2.3----45/.,2");
    samrec.setAttribute(SamUtils.CG_READ_DELTA, "TA");
    samrec.setAttribute(SamUtils.ATTRIBUTE_ALIGNMENT_SCORE, 1);
    samrec.setFlags(67);
    validator.setData(samrec, DnaUtils.encodeString("AGCCCCCNTACGTAAATAAGACATC     ACGATGATCA".replaceAll(" ", "")), FastaUtils.asciiToRawQuality("4316%%68883-56+141663,2.3----45/.,2"));
    validator.setTemplate(DnaUtils.encodeString("AGCCCCCACNCGTAAATAAGACATCTTTTTACGATGATCA"));

    validator.parse();
    assertFalse(validator.getInvalidReason(), validator.isValid());
    assertEquals(validator.getInvalidReason(), "Super cigar alignment score was 3, but AS attribute was 1", validator.getInvalidReason());

    samrec.setAttribute(SamUtils.ATTRIBUTE_ALIGNMENT_SCORE, 3);
    validator.setData(samrec, DnaUtils.encodeString("AGCCCCCNTACGTAAATAAGACATC     ACGATGATCA".replaceAll(" ", "")), FastaUtils.asciiToRawQuality("4316%%68883-56+141663,2.3----45/.,2"));
    validator.setTemplate(DnaUtils.encodeString("AGCCCCCACNCGTAAATAAGACATCTTTTTACGATGATCA"));

    validator.parse();
    assertTrue(validator.getInvalidReason(), validator.isValid());
  }

  public void testOverlapPastStartPosition() throws Exception {

    final SuperCigarValidator validator = new SuperCigarValidator(0);
    final SAMRecord samrec = new SAMRecord(null);

    samrec.setAlignmentStart(3);
    samrec.setCigarString("2=2I16=6N10=");
    samrec.setReadString("ATAAGAAGGAGTGGCACTTCCCTCAGCTCA");
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "2=2I1=4B1X1=1I17=6N10=");
    samrec.setBaseQualityString("20001.1-+,8/0/41373,1751662362");
    samrec.setAttribute(SamUtils.CG_READ_DELTA, "AAGG");
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR_OVERLAP_QUALITY, ").,1/");
    samrec.setAttribute(SamUtils.ATTRIBUTE_ALIGNMENT_SCORE, 6);
    samrec.setFlags(179);

    final String rawread = "TGAGCTGAGG     GAAGTGCCACTCCTTCACTCCTTAT".replaceAll(" ", "");
    final String rawqual = "2632661571     ,37314/0/8,+-1./1,.)10002".replaceAll(" ", "");

    validator.setData(samrec, DnaUtils.encodeString(rawread), FastaUtils.asciiToRawQuality(rawqual));
    validator.setTemplate(DnaUtils.encodeString("TCAT  GAAGGAGTGGCACTTCCACCTGCCTCAGCTCATGCGTGATATCCAGG".replaceAll(" ", "")));
      //                                           ATAAGAAGGAGTGGCACTTCCCTCAGCTCA
      //                                           ATAAG
      //                                         TCA TGAAGGAGTGGCACTTCCACCTGCCTCAGCTCATGCGTGATATCCAGG
      //                                          GAGTGAAGGAGTGGCACTTC      CCTCAGCTCA


    validator.parse();
    assertTrue(validator.getInvalidReason(), validator.isValid());
  }
}
