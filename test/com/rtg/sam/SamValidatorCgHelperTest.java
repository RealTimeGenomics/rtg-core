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

import java.util.Arrays;

import com.rtg.mode.DnaUtils;
import com.rtg.reader.FastaUtils;

import htsjdk.samtools.SAMRecord;

import junit.framework.TestCase;


/**
 */
public class SamValidatorCgHelperTest extends TestCase {

  public void testQual() {
    final SAMRecord samrec = new SAMRecord(null);
    samrec.setAlignmentStart(1);
    samrec.setCigarString("23M7N10M");
    samrec.setReadString("GCTGACCGCCAAAGGTGAGCAACATGAGGTGGC");
    samrec.setBaseQualityString("431%68883-56+141663,2.3-45/.,2553");
    samrec.setAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES, 1);
    samrec.setAttribute(SamUtils.ATTRIBUTE_CG_RAW_READ_INSTRUCTIONS, "3S2G28S");
    samrec.setAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_BASES, "GAGA");
    samrec.setAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_QUALITY, "6%");
    samrec.setReadNegativeStrandFlag(true);
    samrec.setFlags(179);

    byte[] blah = SamValidatorCgHelper.expandCgCigarQualities(samrec.getBaseQualities(), new byte[35], samrec.getStringAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_QUALITY), new int[] {3, 2, 28}, samrec.getFirstOfPairFlag(), samrec.getReadNegativeStrandFlag(), true);
    byte[] exp = FastaUtils.asciiToRawQuality("3552,./54-3.2,366141+65-38886%%6134");
    assertTrue(Arrays.toString(exp) + "\n" + Arrays.toString(blah), Arrays.equals(exp, blah));

    samrec.setCigarString("23M7N10M");
    samrec.setReadString("TTTGTAGGTCGGATAAGGCGTTCATCCGACACG");
    samrec.setBaseQualityString("431%68883-56+141663,2.3-45/.,2553");
    samrec.setAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES, 0);
    samrec.setAttribute(SamUtils.ATTRIBUTE_CG_RAW_READ_INSTRUCTIONS, "3S2G28S");
    samrec.setAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_BASES, "GTGT");
    samrec.setAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_QUALITY, "6%");
    samrec.setReadNegativeStrandFlag(false);
    samrec.setFlags(67);

    blah = SamValidatorCgHelper.expandCgCigarQualities(samrec.getBaseQualities(), new byte[35], samrec.getStringAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_QUALITY), new int[] {3, 2, 28}, samrec.getFirstOfPairFlag(), samrec.getReadNegativeStrandFlag(), true);
    assertTrue(Arrays.equals(FastaUtils.asciiToRawQuality("4316%%68883-56+141663,2.3-45/.,2553"), blah));

    samrec.setCigarString("10M6N24M");
    samrec.setReadString("ACATTCATCACAGACCTGGGCCTGCTGGGCCCCA");
    samrec.setBaseQualityString(".531337242.875&3158.,2+/./4/-,4550");
    samrec.setAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES, 0);
    samrec.setAttribute(SamUtils.ATTRIBUTE_CG_RAW_READ_INSTRUCTIONS, "29S1G4S");
    samrec.setAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_BASES, "CC");
    samrec.setAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_QUALITY, "2");
    samrec.setFlags(131);

    blah = SamValidatorCgHelper.expandCgCigarQualities(samrec.getBaseQualities(), new byte[35], samrec.getStringAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_QUALITY), new int[] {29, 1, 4}, samrec.getFirstOfPairFlag(), samrec.getReadNegativeStrandFlag(), true);
    exp = FastaUtils.asciiToRawQuality(".531337242.875&3158.,2+/./4/-,24550");
    assertTrue(Arrays.toString(exp) + "\n" + Arrays.toString(blah), Arrays.equals(exp, blah));

    samrec.setCigarString("10M5N23M");
    samrec.setReadString("TCAGCACTTCGATAGAGGTTTTCACCACGTTCC");
    samrec.setBaseQualityString("431%68883-56+141663,2.3-45/.,2553");
    samrec.setAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES, 0);
    samrec.setAttribute(SamUtils.ATTRIBUTE_CG_RAW_READ_INSTRUCTIONS, "28S2G3S");
    samrec.setAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_BASES, "GTGT");
    samrec.setAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_QUALITY, "6%");
    samrec.setReadNegativeStrandFlag(true);
    samrec.setFlags(115);

    blah = SamValidatorCgHelper.expandCgCigarQualities(samrec.getBaseQualities(), new byte[35], samrec.getStringAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_QUALITY), new int[] {28, 2, 3}, samrec.getFirstOfPairFlag(), samrec.getReadNegativeStrandFlag(), true);
    assertTrue(Arrays.equals(FastaUtils.asciiToRawQuality("355%62,./54-3.2,366141+65-38886%134"), blah));

  }

  public void testSomeCGStuff() {
    final SAMRecord samrec = new SAMRecord(null);

    samrec.setAlignmentStart(1);
    samrec.setCigarString("20M5N10M");
    samrec.setReadString("GAGGCCGAGGCAGGCGGATCGTCAGGAGTT");
    samrec.setAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES, 0);

    assertTrue(SamValidatorCgHelper.matchesCg(DnaUtils.encodeArray("gaggccgagg     caggcggatcgtcaggagtt".replaceAll(" ", "").getBytes()), null, samrec, false));

    samrec.setAlignmentStart(1);
    samrec.setCigarString("10M5N20M");
    samrec.setReadString("AACTCCTGACGATCCGCCTGCCTCGGCCTC");
    samrec.setAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES, 0);
    samrec.setReadNegativeStrandFlag(true);

    assertTrue(SamValidatorCgHelper.matchesCg(DnaUtils.encodeArray("gaggccgagg     caggcggatcgtcaggagtt".replaceAll(" ", "").getBytes()), null, samrec, false));

    samrec.setCigarString("23M7N10M");
    samrec.setReadString("TTTGTAGGTCGGATAAGGCGTTCATCCGACACG");
    samrec.setBaseQualityString("431%68883-56+141663,2.3-45/.,2553");
    samrec.setAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES, 0);
    samrec.setAttribute(SamUtils.ATTRIBUTE_CG_RAW_READ_INSTRUCTIONS, "3S2G28S");
    samrec.setAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_BASES, "GTGT");
    samrec.setAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_QUALITY, "6%");
    samrec.setReadNegativeStrandFlag(false);
    samrec.setFlags(67);
    assertTrue(SamValidatorCgHelper.matchesCg(DnaUtils.encodeArray("tttgtgtaggtcggataaggcgttc     atccgacacg".replaceAll(" ", "").getBytes()), FastaUtils.asciiToRawQuality("4316%%68883-56+141663,2.3-45/.,2553"), samrec, false));
    //              "431  %68883-56+141663,2.3-45/.,2553"

    // 353     179     NC_000913       4631900 255     23M5N3M1D7M     =       4632182 281     GCTGACCGCCAAAGGTGAGCAACATGAGGTGGC       *       AS:i:2  NM:i:1  MQ:i:255        GS:Z:GAG

    samrec.setCigarString("23M7N10M");
    samrec.setReadString("GCTGACCGCCAAAGGTGAGCAACATGAGGTGGC");
    samrec.setBaseQualityString("431%68883-56+141663,2.3-45/.,2553");
    samrec.setAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES, 1);
    samrec.setAttribute(SamUtils.ATTRIBUTE_CG_RAW_READ_INSTRUCTIONS, "3S2G28S");
    samrec.setAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_BASES, "GAGA");
    samrec.setAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_QUALITY, "6%");
    samrec.setReadNegativeStrandFlag(true);
    samrec.setFlags(179);
    assertTrue(SamValidatorCgHelper.matchesCg(DnaUtils.encodeArray("gccacctcat     gttgctcacctttggcggtctcagc".replaceAll(" ", "").getBytes()), FastaUtils.asciiToRawQuality("3552,./54-3.2,366141+65-38886%%6134"), samrec, false));

    //Read doesn't match expected value from SDF 861  131     NC_000913       4606315 255     10M5N24M        NC_000913       4605888 -427    GGGCGAGCGTCAACTTTCAGTTAACGCAAAACGG     * AS:i:0  NM:i:0  MQ:i:255        GS:Z:AA GC:Z:29S1G4S    XA:i:0  IH:i:1
    samrec.setCigarString("10M5N24M");
    samrec.setReadString("GGGCGAGCGTCAACTTTCAGTTAACGCAAAACGG");
    samrec.setAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES, 0);
    samrec.setAttribute(SamUtils.ATTRIBUTE_CG_RAW_READ_INSTRUCTIONS, "29S1G4S");
    samrec.setAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_BASES, "AT");
    samrec.setAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_QUALITY, "6");
    samrec.setReadNegativeStrandFlag(false);
    samrec.setFlags(131);
    assertTrue(SamValidatorCgHelper.matchesCg(DnaUtils.encodeArray("gggcgagcgt     caactttcagttaacgcaaatacgg".replaceAll(" ", "").getBytes()), null, samrec, false));

    //Read doesn't match expected value from SDF 572  115     NC_000913       855891  255     10M5N16M1D8M    NC_000913       855206  -683    AACTACGATGATGACGGCGGCGACCCGGCGCGTG     *    AS:i:2  NM:i:1  MQ:i:255        GS:Z:GA GC:Z:29S1G4S    XA:i:4  IH:i:1

    samrec.setCigarString("10M5N16M1D8M");
    samrec.setReadString("AACTACGATGATGACGGCGGCGACCCGGCGCGTG");
    samrec.setAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES, 1);
    samrec.setAttribute(SamUtils.ATTRIBUTE_CG_RAW_READ_INSTRUCTIONS, "29S1G4S");
    samrec.setAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_BASES, "GA");
    samrec.setReadNegativeStrandFlag(true);
    samrec.setFlags(115);
    assertTrue(SamValidatorCgHelper.matchesCg(DnaUtils.encodeArray("cacgtcgccgggtcgccgccgtcat     catcgtagtt".replaceAll(" ", "").getBytes()), null, samrec, false));


    //Read doesn't match expected value from SDF 703  115     NC_000913       4584340 255     10M5N23M        NC_000913       4584024 -314    TCAGCACTTCGATAGAGGTTTTCACCACGTTCC      *    AS:i:0  NM:i:0  MQ:i:255        GS:Z:GTGT       GC:Z:28S2G3S    XA:i:0  IH:i:1
    samrec.setCigarString("10M5N23M");
    samrec.setReadString("TCAGCACTTCGATAGAGGTTTTCACCACGTTCC");
    samrec.setAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES, 0);
    samrec.setAttribute(SamUtils.ATTRIBUTE_CG_RAW_READ_INSTRUCTIONS, "28S2G3S");
    samrec.setAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_BASES, "GTGT");
    samrec.setReadNegativeStrandFlag(true);
    samrec.setFlags(115);
    assertTrue(SamValidatorCgHelper.matchesCg(DnaUtils.encodeArray("ggaacacgtggtgaaaacctctatc     gaagtgctga".replaceAll(" ", "").getBytes()), null, samrec, false));


    //15404   131     paolo-bac       100219  255     10M6N24M        =       99888   -331    ACATTCATCACAGACCTGGGCCTGCTGGGCCCCA      .531337242.875&3158.,2+/./4/-,4550      AS:i:0    NM:i:0  MQ:i:255        GS:Z:CC GC:Z:29S1G4S    GQ:Z:2  XA:i:2  IH:i:1

    samrec.setCigarString("10M6N24M");
    samrec.setReadString("ACATTCATCACAGACCTGGGCCTGCTGGGCCCCA");
    samrec.setBaseQualityString(".531337242.875&3158.,2+/./4/-,4550");
    samrec.setAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES, 0);
    samrec.setAttribute(SamUtils.ATTRIBUTE_CG_RAW_READ_INSTRUCTIONS, "29S1G4S");
    samrec.setAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_BASES, "CC");
    samrec.setAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_QUALITY, "2");
    samrec.setFlags(131);
    assertTrue(SamValidatorCgHelper.matchesCg(DnaUtils.encodeArray("acattcatca     cagacctgggcctgctgggccccca".replaceAll(" ", "").getBytes()), FastaUtils.asciiToRawQuality(".531337242.875&3158.,2+/./4/-,24550"), samrec, false));


    samrec.setCigarString("25M6N10M");
    samrec.setReadString("TATTAGGATTGAGACTGGTAAAATGGNCCACCAAG");
    samrec.setBaseQualityString("/725361840-525251.68,0,.52!*/*/54/2");
    samrec.setAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES, 2);
    samrec.setAttribute(SamUtils.ATTRIBUTE_CG_RAW_READ_INSTRUCTIONS, null);
    samrec.setAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_BASES, null);
    samrec.setAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_QUALITY, null);
    samrec.setFlags(67);
    assertTrue(SamValidatorCgHelper.matchesCg(DnaUtils.encodeArray("tattaggatt     gagactggtaaaatggnccaccaag".replaceAll(" ", "").getBytes()), FastaUtils.asciiToRawQuality("/725361840-525251.68,0,.52!*/*/54/2"), samrec, false));

    //Read doesn't match expected value from SDF 40928        131     paolo-bac       411     255     10M5N16M1I3M2I1M        paolo-bac       93      -318    TTCTCCATTCGAGACGTTGTGAAT    GTGGACTTG       .333,0/60..855041063+.,1014/0*#45       AS:i:7  NM:i:6  MQ:i:255        GS:Z:ACAC       GC:Z:28S2G3S    GQ:Z:06 XA:i:8  IH:i:1

    samrec.setCigarString("10M5N16M1I3M2I1M");
    samrec.setReadString("TTCTCCATTCGAGACGTTGTGAATGTGGACTTG");
    samrec.setBaseQualityString(".333,0/60..855041063+.,1014/0*#45");
    samrec.setAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES, 6);
    samrec.setAttribute(SamUtils.ATTRIBUTE_CG_RAW_READ_INSTRUCTIONS, "28S2G3S");
    samrec.setAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_BASES, "ACAC");
    samrec.setAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_QUALITY, "06");
    samrec.setFlags(131);
    assertTrue(SamValidatorCgHelper.matchesCg(DnaUtils.encodeArray("ttctccattc     gagacgttgtgaatgtggacacttg".replaceAll(" ", "").getBytes()), FastaUtils.asciiToRawQuality(".333,0/60..855041063+.,1014/0*06#45"), samrec, false));
  }

  public void testSomeBadCGStuff() {
    final SAMRecord samrec = new SAMRecord(null);

    samrec.setCigarString("10=5N10=2I2=1X1=1X7=");
    samrec.setReadString("TAAAATATTTTAGTCTTTTACCTAAATCTTCAAA");
    samrec.setBaseQualityString("7/8998-9:7;8<;9<;<;96777667217:887");
    samrec.setAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES, 5);
    samrec.setAttribute(SamUtils.ATTRIBUTE_CG_RAW_READ_INSTRUCTIONS, "30S1G3S");
    samrec.setAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_BASES, "CT");
    samrec.setAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_QUALITY, ":");
    samrec.setReadNegativeStrandFlag(true);
    samrec.setFirstOfPairFlag(true);
    samrec.setFlags(115);
    assertFalse(SamValidatorCgHelper.matchesCg(DnaUtils.encodeArray("tttagaagatttaggtaaaagacta     aaatatttta".replaceAll(" ", "").getBytes()), FastaUtils.asciiToRawQuality("788::71276677769;<;<9;<8;7:9-8998/7"), samrec, false));
    assertTrue(SamValidatorCgHelper.matchesCg(DnaUtils.encodeArray("tttagaagatttaggtaaaagacta     aaatatttta".replaceAll(" ", "").getBytes()), FastaUtils.asciiToRawQuality("788::71276677769;<;<9;<8;7:9-8998/7"), samrec, true));

    samrec.setCigarString("8=1I1X2=1X11=6N10=");
    samrec.setReadString("ATATGTATATATCAGAGGTTAATAGTTATGAGGA");
    samrec.setBaseQualityString("295952025637769;;<;:8.:78888988884");
    samrec.setAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES, 4);
    samrec.setAttribute(SamUtils.ATTRIBUTE_CG_RAW_READ_INSTRUCTIONS, "3S1G30S");
    samrec.setAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_BASES, "GT");
    samrec.setAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_QUALITY, "6");
    samrec.setReadNegativeStrandFlag(false);
    samrec.setFirstOfPairFlag(true);
    samrec.setFlags(67);
    assertFalse(SamValidatorCgHelper.matchesCg(DnaUtils.encodeArray("ATAGTGTATATATCAGAGGTTAATA     GTTATGAGGA".replaceAll(" ", "").getBytes()), FastaUtils.asciiToRawQuality("2956952025637769;;<;:8.:78888988884"), samrec, false));
    assertTrue(SamValidatorCgHelper.matchesCg(DnaUtils.encodeArray("ATAGTGTATATATCAGAGGTTAATA     GTTATGAGGA".replaceAll(" ", "").getBytes()), FastaUtils.asciiToRawQuality("2956952025637769;;<;:8.:78888988884"), samrec, true));

  }
}
