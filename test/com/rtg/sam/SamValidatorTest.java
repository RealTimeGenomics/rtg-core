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

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;

import com.rtg.mode.DNA;
import com.rtg.mode.DnaUtils;
import com.rtg.ngs.NgsParams;
import com.rtg.ngs.NgsParamsBuilder;
import com.rtg.reader.PrereadArm;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.ReaderUtils;
import com.rtg.reader.SdfId;
import com.rtg.sam.SamValidator.SamStatsVariables;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;

import htsjdk.samtools.SAMRecord;
import junit.framework.TestCase;

/**
 * Tests corresponding class
 */
public class SamValidatorTest extends TestCase {

  public void testSamStatsVariables() {
    final SamStatsVariables total = new SamStatsVariables();
    assertEquals(0, total.mMaxInsertSize);
    assertEquals(Integer.MAX_VALUE, total.mMinInsertSize);
    assertEquals(0, total.mMaxReadLength);
    assertEquals(Integer.MAX_VALUE, total.mMinReadLength);
    assertEquals(0, total.mTotalRecords);
    assertEquals(0, total.mTotalMatches);
    assertEquals(0, total.mTotalMismatches);
    assertEquals(0, total.mUnmappedRecords);
    assertEquals(4, total.mPairOrientations.length);
    assertEquals(0, total.mConsensus);
    assertEquals(0, total.mTotalNt);
    assertEquals(0, total.mTotalLength);
    assertEquals(0, total.mAlignmentScores.size());
    assertEquals(0, total.mInsertSizes.size());
    assertEquals(0, total.mNH.size());

    final SamStatsVariables subunit = new SamStatsVariables();
    subunit.mMaxInsertSize = 3;
    subunit.mMinInsertSize = 2;
    subunit.mMaxReadLength = 70;
    subunit.mMinReadLength = 20;
    subunit.mTotalRecords = 4;
    subunit.mTotalMatches = 5;
    subunit.mTotalMismatches = 6;
    subunit.mUnmappedRecords = 1;
    for (int i = 0; i < subunit.mPairOrientations.length; i++) {
      subunit.mPairOrientations[i] = (i + 1) * 7;
    }
    subunit.mConsensus = 8;
    subunit.mTotalNt = 9;
    subunit.mTotalLength = 10;
    subunit.mAlignmentScores.put(11, 11);
    subunit.mInsertSizes.put(12, 12);
    subunit.mNH.put(13, 13);
    total.addToTotal(subunit);
    total.addToTotal(subunit);
    assertEquals(3, total.mMaxInsertSize);
    assertEquals(2, total.mMinInsertSize);
    assertEquals(70, total.mMaxReadLength);
    assertEquals(20, total.mMinReadLength);
    assertEquals(8, total.mTotalRecords);
    assertEquals(10, total.mTotalMatches);
    assertEquals(12, total.mTotalMismatches);
    assertEquals(2, total.mUnmappedRecords);
    assertEquals(4, total.mPairOrientations.length);
    for (int i = 0; i < total.mPairOrientations.length; i++) {
      assertEquals((i + 1) * 7 * 2, total.mPairOrientations[i]);
    }
    assertEquals(16, total.mConsensus);
    assertEquals(18, total.mTotalNt);
    assertEquals(20, total.mTotalLength);
    assertEquals(1, total.mAlignmentScores.size());
    assertEquals((Integer) 22, total.mAlignmentScores.get(11));
    assertEquals(1, total.mInsertSizes.size());
    assertEquals((Integer) 24, total.mInsertSizes.get(12));
    assertEquals(1, total.mNH.size());
    assertEquals((Integer) 26, total.mNH.get(13));
  }

  public void testPrint2() {

    final int[] countStats = {1, 4, 2, 4, 5, 0, 2};
    final boolean[] paired = {false, true, true, false, false, false, false};

    final ByteArrayOutputStream os = new ByteArrayOutputStream();
    try (PrintStream newOut = new PrintStream(os)) {
      SamValidator.printStats2(countStats, paired, newOut, false);
      newOut.flush();
      final String out = os.toString();

      assertTrue(out.contains("Reads seen: 6 of 7"));
      assertTrue(out.contains("Max records for a read: 5"));
      assertTrue(out.contains("Reads mapped at a single location: 1"));
      assertTrue(out.contains("0 : 1"));
      assertTrue(out.contains("1 : 1"));
      assertTrue(out.contains("2 : 2"));
      assertTrue(out.contains("4 : 2"));
      assertTrue(out.contains("5 : 1"));
    } finally {
      try {
        os.close();
      } catch (final IOException ioe) {
      }
    }
  }

  private byte[] stringToDna(final String s) {
    final byte[] b = new byte[s.length()];
    for (int i = 0; i < s.length(); i++) {
      b[i] = (byte) DNA.valueOf(String.valueOf(s.charAt(i))).ordinal();
    }
    return b;
  }


  private static final String SAM_EXAMPLE1 = ""
    + "@HD\tVN:1.0\tSO:coordinate" + StringUtils.LS
    + "@SQ\tSN:simulatedSequence1\tLN:148" + StringUtils.LS
    + "@SQ\tSN:simulatedSequence2\tLN:134" + StringUtils.LS
    + "@PG\tID:slim VN:v2.0-EAP1.3 build 23142 (2009-12-16)" + StringUtils.LS
    + "3\t83\tsimulatedSequence1\t3\t255\t2M\t=\t15\t0\tAA\t*\tAS:i:0\tNM:i:0\tMQ:i:255\tIH:i:1" + StringUtils.LS
    + "3\t163\tsimulatedSequence1\t15\t255\t2M\t=\t3\t0\tTT\t*\tAS:i:0\tNM:i:0\tMQ:i:255\tIH:i:1" + StringUtils.LS
    + "0\t163\tsimulatedSequence1\t24\t255\t30M\t=\t108\t114\tCGTAGTGGAATCGATGCTAATGAGGACCAG\t*\tAS:i:0\tNM:i:0\tMQ:i:255\tIH:i:1" + StringUtils.LS
    + "2\t163\tsimulatedSequence1\t49\t255\t30M\t=\t102\t83\tACCAGGGTACGAAAAGACACGATAAACTAC\t*\tAS:i:0\tNM:i:0\tMQ:i:255\tIH:i:1" + StringUtils.LS
    + "1\t99\tsimulatedSequence1\t86\t255\t30M\t=\t89\t33\tCGAGAATGGTTACCCCTCCCAGATGGCTCG\t*\tAS:i:0\tNM:i:0\tMQ:i:255\tIH:i:1" + StringUtils.LS
    + "1\t147\tsimulatedSequence1\t89\t255\t30M\t=\t86\t-33\tGAATGGTTACCCCTCCCAGATGGCTCGCCA\t*\tAS:i:0\tNM:i:0\tMQ:i:255\tIH:i:1" + StringUtils.LS
    + "2\t83\tsimulatedSequence1\t102\t255\t30M\t=\t49\t-83\tTCCCAGATGGCTCGCCACACGGCTAACAAG\t*\tAS:i:0\tNM:i:0\tMQ:i:255\tIH:i:1" + StringUtils.LS
    + "0\t83\tsimulatedSequence1\t108\t255\t30M\t=\t24\t-114\tATGGCTCTCCACACGGCTAACAAGTGAGTG\t*\tAS:i:3\tNM:i:3\tMQ:i:255\tIH:i:1" + StringUtils.LS;


  private static final String READS_LEFT = ">read0:simulatedSequence1:138R:S1:I0:D1" + StringUtils.LS
  + "CACTCACTTGTTAGCCGTGTGGAGAGCCAT" + StringUtils.LS
  + ">read1:simulatedSequence1:86:S0:I0:D0" + StringUtils.LS
  + "CGAGAATGGTTACCCCTCCCAGATGGCTCG" + StringUtils.LS
  + ">read2:simulatedSequence1:131R:S0:I0:D0" + StringUtils.LS
  + "CTTGTTAGCCGTGTGGCGAGCCATCTGGGA" + StringUtils.LS
  + ">read3:simulatedSequence1:49:S0:I0:D0" + StringUtils.LS
  + "TT" + StringUtils.LS;


  private static final String READS_RIGHT = ">read0:simulatedSequence1:24:S0:I0:D0" + StringUtils.LS
  + "CGTAGTGGAATCGATGCTAATGAGGACCAG" + StringUtils.LS
  + ">read1:simulatedSequence1:118R:S0:I0:D0" + StringUtils.LS
  + "TGGCGAGCCATCTGGGAGGGGTAACCATTC" + StringUtils.LS
  + ">read2:simulatedSequence1:49:S0:I0:D0" + StringUtils.LS
  + "ACCAGGGTACGAAAAGACACGATAAACTAC" + StringUtils.LS
  + ">read3:simulatedSequence1:49:S0:I0:D0" + StringUtils.LS
  + "TT" + StringUtils.LS;

  private static final String TEMPLATE = ">simulatedSequence1" + StringUtils.LS
  + "CCAACGTGAATATGTTTCCTCACCGTAGTGGAATCGATGCTAATGAGGACCAGGGTACGAAAAGACACGATAAACTACCTAGTCACGAGAATGGTTACCCCTCCCAGATGGCTCGCCACACGGCTAACAAGTGAGGTGGCACCGTAAG" + StringUtils.LS
  + ">simulatedSequence2" + StringUtils.LS
  + "TTTGTAACGCAGGCATTCGTGCGATAAGGTTGAGGGCCAAGGACGGACACGGCCGATCCGGACATTGGCGCTCTGTTTACCCTCTAAGTCCGTTAGGGGGGATTAGATGTCTGACTATCCCGCAAAGTGGGTTG" + StringUtils.LS;


  public void testWhole() throws Exception {
    Diagnostic.setLogStream();
    try (final TestDirectory tmpDir = new TestDirectory("samval")) {
      final ByteArrayOutputStream baos = new ByteArrayOutputStream();
      try (PrintStream ps = new PrintStream(baos)) {

        final File samFile = new File(tmpDir, "sam.sam");
        FileUtils.stringToFile(SAM_EXAMPLE1, samFile);
        final Collection<File> files = new ArrayList<>();
        files.add(samFile);
        final SamValidator sv = new SamValidator(ps, ps, true, true, true, true, true, getParams(false), true);

        final File templateDir = new File(tmpDir, "templ");
        ReaderTestUtils.getReaderDNA(TEMPLATE, templateDir, null);
        ReaderTestUtils.createPairedReaderDNA(READS_LEFT, READS_RIGHT, tmpDir, null);
        sv.checkSAMAlign(templateDir, files, ReaderUtils.getLeftEnd(tmpDir), ReaderUtils.getRightEnd(tmpDir));
        ps.flush();

        final String[] expected = {
                                          "Stats for file: " + samFile.getPath(),
                                          "Total records: 8",
                                          "Unmapped records: 0",
                                          "Coverage: 1.2432",
                                          "Consensus: 0.9946",
                                          "% accuracy compared to reference: 98.40",
                                          "FF: 0",
                                          "RF: 1",
                                          "FR: 3",
                                          "RR: 0",
                                          "Min/Max insert size: 0/114",
                                          "Min/Max read length: 2/30",
                                          "Reads seen: 4 of 4",
                                          "Reads mapped at a single location: 4",
                                          "Max records for a read: 2",
                                          "Record counts by read:",
                                          "2 : 4",
                                          "Distribution of alignment score:",
                                          "#tag\tscore\tcount",
                                          "AS:\t0\t7",
                                          "AS:\t3\t1",
                                          "Distribution of read hits:",
                                          "#tag\thits\tcount",
                                          "IH:\t1\t8",
                                          "Distribution of insert size:",
                                          "#tag\tsize\tcount",
                                          "IS:\t0\t1",
                                          "IS:\t33\t1",
                                          "IS:\t83\t1",
                                          "IS:\t114\t1"
        };


        assertFalse(baos.toString().contains("Read doesn't match expected value from SDF"));
        assertFalse(baos.toString().contains("Alignment score disparity"));

        TestUtils.containsAll(baos.toString(), expected);
      }
    }
  }

  public void testBadInputs() throws Exception {
    try (final TestDirectory tmpDir = new TestDirectory("samval")) {
      final ByteArrayOutputStream baos = new ByteArrayOutputStream();
      final PrintStream ps = new PrintStream(baos);
      Diagnostic.setLogStream(ps);
      try {
        final File samFile = new File(tmpDir, "sam.sam");
        FileUtils.stringToFile(SAM_EXAMPLE1, samFile);
        final Collection<File> files = new ArrayList<>();
        files.add(samFile);
        final SamValidator sv = new SamValidator(null, null, true, true, false, false, false, getParams(false), false);

        final File templateDir = new File(tmpDir, "templ");
        ReaderTestUtils.getReaderDNA(TEMPLATE, templateDir, null);
        final File leftReadsDir = new File(tmpDir, "left");
        ReaderTestUtils.getReaderDNA(READS_LEFT, leftReadsDir, true, new SdfId(27));
        final File rightReadsDir = new File(tmpDir, "right");
        ReaderTestUtils.getReaderDNA(READS_RIGHT, rightReadsDir, true, new SdfId(42));

        try {
          sv.checkSAMAlign(templateDir, files, leftReadsDir, rightReadsDir);
          fail();
        } catch (final NoTalkbackSlimException ntse) {
          ntse.logException();
          ps.flush();
          assertTrue(baos.toString().contains("Error: Problem reading file: \"Left and right reads have different GUIDs - are from different runs.\""));
        }

        final File leftCGReadsDir = new File(tmpDir, "leftcg");
        ReaderTestUtils.getReaderDNAFastqCG("", leftCGReadsDir, PrereadArm.LEFT);
        try {
          System.out.println(Arrays.toString(leftCGReadsDir.listFiles()));
          sv.checkSAMAlign(templateDir, files, leftCGReadsDir, rightReadsDir);
          fail();
        } catch (final NoTalkbackSlimException ntse) {
          ntse.logException();
          ps.flush();
          assertTrue(baos.toString().contains("Error: Problem reading file: \"Left reads are CG data, but right reads are not.\""));
        }

        final File sam = new File(tmpDir, "samvalnosort.sam");
        final String tab = "\t";
        final String samHeader = ""
          + "@HD" + tab + "VN:1.0" + tab + "SO:unsorted" + StringUtils.LS
          + "@SQ" + tab + "SN:gi0" + tab + "LN:30" + StringUtils.LS
          + "@SQ" + tab + "SN:gi1" + tab + "LN:30" + StringUtils.LS
          + "aa" + tab + "0" + tab + "gi0" + tab + "2" + tab + "255" + tab + "10M" + tab + "*" + tab + "0" + tab + "0" + tab + "AAAAAAAAAA" + tab +  "IB7?*III<I" + tab + "AS:i:0" + tab + "IH:i:1" + StringUtils.LS
          ;
        FileUtils.stringToFile(samHeader, sam);
        files.clear();
        files.add(sam);
        final File rightReadsDir2 = new File(tmpDir, "right2");
        ReaderTestUtils.getReaderDNA(READS_RIGHT, rightReadsDir2, true, new SdfId(27));
        try {
          sv.checkSAMAlign(templateDir, files, leftReadsDir, rightReadsDir2);
          fail();
        } catch (final NoTalkbackSlimException ntse) {
          ntse.logException();
          ps.flush();
          assertTrue(baos.toString().contains("Error: Problem reading file: \"SAM file must be sorted.\""));
        }

      } finally {
        Diagnostic.setLogStream();
        ps.close();
      }
    }
  }


  public void testAlignment() {
    String template = "ATGTTCTGTGCATCTTCCCTTACCTGGGGCCCTCACTGAGTGGGTCCTCCATGGGTGACTGGTGA";

    final SAMRecord samrec = new SAMRecord(null);
    try (PrintStream out = new PrintStream(new ByteArrayOutputStream())) {
      samrec.setAlignmentStart(6);
      samrec.setCigarString("4M1D6M1D15M");
      samrec.setReadString("ATGTCATCTTCCTCTTCTGGGGCNT");
      samrec.setBaseQualityString("55663,0/5534.2,898256678.3887755562");
      samrec.setAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES, 8);
      final SamValidator sv = new SamValidator(out, out, false, false, false, false, false, getParams(false), true);
      assertEquals(9, sv.isAtExpectedRef(stringToDna(template), samrec, new PileUp(template.length())));

      samrec.setBaseQualityString("55663,0/5534.2,898256678.");

      assertFalse(sv.matchesRawRead(DnaUtils.encodeArray("atgtcatcttcctcttctggggcnt".getBytes()), DnaUtils.fastqToPhred("55663,0/5534.2,898256678.3887755562"), samrec, false, false));
      assertTrue(sv.matchesRawRead(DnaUtils.encodeArray("atgtcatcttcctcttctggggcnt".getBytes()), DnaUtils.fastqToPhred("55663,0/5534.2,898256678."), samrec, false, false));

      template = "GGCACACCTGGAGCCAGCCACCCGCTGGGTCGCACATGGATCTGGTGATATTATTGATAAT";
      samrec.setReadString("GGCACACCTGCCACCCGCTGGGTCGCACATGGA");
      samrec.setAlignmentStart(1);
      samrec.setCigarString("10M7N23M");
      samrec.setAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES, 0);
      samrec.setAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_BASES, "ATAT");
      samrec.setAttribute(SamUtils.ATTRIBUTE_CG_RAW_READ_INSTRUCTIONS, "28S2G3S");  //28 = 35 - 4 - 3
      samrec.setFlags(115);
      assertEquals(0, sv.isAtExpectedRef(stringToDna(template), samrec, new PileUp(template.length())));

      template = "AACATGACTAAAGTACGTAATTGCGTTCTTGATGCACTTTC";
      samrec.setReadString("AACATGACTAAAGTACGTAATTGCTGTGCACTTT");
      samrec.setAlignmentStart(1);
      samrec.setCigarString("24M5N2M1D8M");
      samrec.setAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES, 1);
      samrec.setAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_BASES, "TT");
      samrec.setAttribute(SamUtils.ATTRIBUTE_CG_RAW_READ_INSTRUCTIONS, "4S1G29S");  //29 = 35 - 4 - 2
      samrec.setFlags(179);
      assertEquals(2, sv.isAtExpectedRef(stringToDna(template), samrec, new PileUp(template.length())));

      template = "GCTGGCGTCTGGCTGGGTCGTTGAAACCGCAGGGGACATCT";
      samrec.setReadString("GCTGGCGTCTGGCTGGGTCGTTGAGGGGACATC");
      samrec.setAlignmentStart(1);
      samrec.setCigarString("23M7N10M");
      samrec.setAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES, 0);
      samrec.setAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_BASES, "GGGG");
      samrec.setAttribute(SamUtils.ATTRIBUTE_CG_RAW_READ_INSTRUCTIONS, "3S2G28S");
      samrec.setFlags(67);
      assertEquals(0, sv.isAtExpectedRef(stringToDna(template), samrec, new PileUp(template.length())));

      template = "GGGCGAGCGTTATCTCAACTTTCAGTTAACGCAAAACGGCAAAATGGCGGC";
      samrec.setReadString("GGGCGAGCGTCAACTTTCAGTTAACGCAAAACGG");
      samrec.setAlignmentStart(1);
      samrec.setCigarString("10M5N24M");
      samrec.setAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES, 0);
      samrec.setAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_BASES, "AA");
      samrec.setAttribute(SamUtils.ATTRIBUTE_CG_RAW_READ_INSTRUCTIONS, "29S1G4S");
      samrec.setFlags(131);
      assertEquals(0, sv.isAtExpectedRef(stringToDna(template), samrec, new PileUp(template.length())));

      template = "CGACGGCGGTGGGATTGCTTCACTATGGGAAAGAGTCACATCTTAACGGTGAAGCTGAAGT";
      samrec.setReadString("CGACGGCGGTTGCTTCACTATGGGAAAGAGTCCA");
      samrec.setAlignmentStart(1);
      samrec.setCigarString("10M5N22M1D2M");
      samrec.setAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES, 1);
      samrec.setAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_BASES, "GG");
      samrec.setAttribute(SamUtils.ATTRIBUTE_CG_RAW_READ_INSTRUCTIONS, "29S1G4S");
      samrec.setFlags(131);
      assertEquals(2, sv.isAtExpectedRef(stringToDna(template), samrec, new PileUp(template.length())));

      final PileUp p = new PileUp(template.length());
      template = "ACACTGAGAAAGGGTATCACTTCGTTAGGGTGGCGGCAGC";
      samrec.setReadString("ACCTGAGAAAGGGTATCACTTCGTTGGCGGCAGC");
      samrec.setAlignmentStart(1);
      samrec.setCigarString("2M1D22M5N10M");
      samrec.setAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES, 1);
      samrec.setAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_BASES, "GG");
      samrec.setAttribute(SamUtils.ATTRIBUTE_CG_RAW_READ_INSTRUCTIONS, "4S1G29S");
      samrec.setFlags(67);
      assertEquals(2, sv.isAtExpectedRef(stringToDna(template), samrec, p));

      template = "CTGCTGCACACAACCCCAGTAAATATCGTCGAGGGCCGCC";
      samrec.setReadString("CTGCTGCACAAGTAAATATCGTCGAGGGCTGCC");
      samrec.setAlignmentStart(1);
      samrec.setCigarString("10M7N23M");
      samrec.setAttribute(SamUtils.ATTRIBUTE_ALIGNMENT_SCORE, 1);
      samrec.setAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES, 1);
      samrec.setAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_BASES, "CC");
      samrec.setAttribute(SamUtils.ATTRIBUTE_CG_RAW_READ_INSTRUCTIONS, "29S1G4S");
      samrec.setFlags(131);
      assertEquals(1, sv.isAtExpectedRef(stringToDna(template), samrec, p));

      assertEquals(47, p.consensus());
      assertEquals(67, p.total());
      assertEquals(1.098360655737705, p.coverage());
    }
  }

  public void testCGSoftClip() {

    final String template = "AGCCCACACGTTCCCCTTAAATAAGACATCACGATG";

    final SAMRecord samrec = new SAMRecord(null);

    samrec.setAlignmentStart(1);
    samrec.setCigarString("10=7N19=5S");

    //AGCCCACACGTAAATAAGACATCACGATGGATCA
    samrec.setReadString("AGCCCACACGTAAATAAGACATCACGATGGATCA");
    samrec.setBaseQualityString(":9999:89;:*;<:2==:<;;::;;68686;8;:");
    samrec.setAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES, 0);
    samrec.setAttribute(SamUtils.ATTRIBUTE_CG_RAW_READ_INSTRUCTIONS, "29S1G4S");
    samrec.setAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_BASES, "GG");
    samrec.setAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_QUALITY, ":");

    try (PrintStream out = new PrintStream(new ByteArrayOutputStream())) {
      final SamValidator sv = new SamValidator(out, out, false, false, false, false, false, getParams(false), false);
      final byte[] templ = stringToDna(template);
      final PileUp p = new PileUp(template.length());
      assertEquals(Integer.MIN_VALUE, sv.isAtExpectedRef(templ, samrec, p));
    }
  }

  public void testNewMatchMismatchStuff() {
    final String template = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG";

    final SAMRecord samrec = new SAMRecord(null);

    samrec.setAlignmentStart(6);
    samrec.setCigarString("35M");
    samrec.setReadString("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG");
    samrec.setBaseQualityString("55663,0/5534.2,898256678.3887755562");
    samrec.setAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES, 0);
    try (PrintStream out = new PrintStream(new ByteArrayOutputStream())) {
      final SamValidator sv = new SamValidator(out, out, false, false, false, false, false, getParams(false), true);
      final byte[] templ = stringToDna(template);
      final PileUp p = new PileUp(template.length());
      assertEquals(0, sv.isAtExpectedRef(templ, samrec, p));
      samrec.setCigarString("35=");
      assertEquals(35, p.consensus());
      assertEquals(0.8333333333333334, p.coverage());
      assertEquals(35, p.total());
      assertEquals(0, sv.isAtExpectedRef(templ, samrec, p));
      assertEquals(70, p.consensus());
      assertEquals(1.6666666666666667, p.coverage());
      assertEquals(70, p.total());

      samrec.setReadString("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
      samrec.setCigarString("35X");
      assertEquals(-1, sv.isAtExpectedRef(templ, samrec, p));

      assertEquals(70, p.consensus());
      assertEquals(2.5, p.coverage());
      assertEquals(105, p.total());
    }
  }

  public void testNewMatchMismatchStuffNs() {
    try (MemoryPrintStream mps = new MemoryPrintStream()) {
      final SAMRecord samrec = new SAMRecord(null);
      final String template = "GGGGGGNNGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG";
      samrec.setReadString("        NGNGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG".replaceAll(" ", ""));
      samrec.setAlignmentStart(6);
      samrec.setCigarString("35M");
      samrec.setBaseQualityString("55663,0/5534.2,898256678.3887755562");
      samrec.setAttribute(SamUtils.ATTRIBUTE_ALIGNMENT_SCORE, 0);
      SamValidator sv = new SamValidator(mps.printStream(), mps.printStream(), false, false, false, false, false, getParams(false), true);
      final byte[] templ = stringToDna(template);
      final PileUp p = new PileUp(template.length());
      assertEquals(mps.toString(), 0, sv.isAtExpectedRef(templ, samrec, p));   //test old cigar with ns as matches

      samrec.setCigarString("3X32=");
      assertEquals(0, sv.isAtExpectedRef(templ, samrec, p));   //test new cigar with ns as matches

      mps.reset();
      samrec.setCigarString("35=");
      assertEquals(0, sv.isAtExpectedRef(templ, samrec, p)); //Ns can be mismatches OR matches now.

      sv = new SamValidator(mps.printStream(), mps.printStream(), false, false, false, false, false, getParams(true), false);
      samrec.setAttribute(SamUtils.ATTRIBUTE_ALIGNMENT_SCORE, 3);
      samrec.setCigarString("35M");
      assertEquals(3, sv.isAtExpectedRef(templ, samrec, p));   //test old cigar with ns as mismatches

      samrec.setCigarString("3X32=");
      assertEquals(3, sv.isAtExpectedRef(templ, samrec, p));   //test new cigar with ns as mismatches

      mps.reset();
      samrec.setCigarString("35=");
      assertEquals(3, sv.isAtExpectedRef(templ, samrec, p)); //Ns can be mismatches OR matches now.
    }
  }

  public void testCgGotohProb() {

    final SAMRecord samrec = new SAMRecord(null);
    samrec.setAlignmentStart(3);
    samrec.setCigarString("2=2I16=6N10=");
    samrec.setReadString("     ATAAGAAGGAGTGGCACTTCC      CTCAGCTCA".replaceAll(" ", ""));
    final String template = "TCAT  GAAGGAGTGGCACTTCCACCTGCCTCAGCTCATGCGTGATATCCAGG".replaceAll(" ", "");
      // orig read             ATAAG
    //                        GA
    //                         GT  GAAGGAGTGGCACTTC
    samrec.setBaseQualityString("20001.1-+,8/0/41373,1751662362");
    samrec.setAttribute(SamUtils.ATTRIBUTE_ALIGNMENT_SCORE, 6);
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "2=2I1=4B1X1=1I17=6N10=");
    samrec.setAttribute(SamUtils.CG_OVERLAP_QUALITY, ").,1/");
    samrec.setAttribute(SamUtils.CG_READ_DELTA, "AAGG");
    samrec.setFlags(179);

//    final String rawread = "TGAGCTGAGG.....GAAGTGCCACTCCTTCACTCCTTAT";
//    final String rawqual = "2632661571     ,37314/0/8,+-1./1,.)10002".replaceAll(" ", "");
//    Expected match, SDF read=A, template=T, 9572    179     paolo-bac       9150    55      2=2I16=6N10=    paolo-bac   9476    329     ATAAGAAGGAGTGGCACTTCCCTCAGCTCA  20001.1-+,8/0/41373,1751662362  AS:i:6  MQ:i:255   XU:Z:2=2I1=4B1X1=1I17=6N10=      XQ:Z:).,1/      XR:Z:AAGG       XA:i:7  IH:i:1  NH:i:1
//    Read doesn't match expected value from SDF 9572 179     paolo-bac       9150    55      2=2I16=6N10=    paolo-bac   9476    329     ATAAGAAGGAGTGGCACTTCCCTCAGCTCA  20001.1-+,8/0/41373,1751662362  AS:i:6  MQ:i:255   XU:Z:2=2I1=4B1X1=1I17=6N10=      XQ:Z:).,1/      XR:Z:AAGG       XA:i:7  IH:i:1  NH:i:1

    try (MemoryPrintStream mps = new MemoryPrintStream()) {
      final SamValidator sv = new SamValidator(mps.printStream(), mps.printStream(), true, false, false, false, false, getParams(false), false);
      final byte[] templ = stringToDna(template);
      final PileUp p = new PileUp(template.length());

      assertEquals(mps.toString(), 3, sv.isAtExpectedRef(templ, samrec, p));  //note this is the CIGAR score, not the super cigar score
    }
  }

  private NgsParams getParams(boolean nsAsMismatches) {
    return new NgsParamsBuilder().gapOpenPenalty(1)
        .gapExtendPenalty(1)
        .substitutionPenalty(1)
        .unknownsPenalty(nsAsMismatches ? 1 : 0).create();
  }

  public void testCgGotohProb2() {

//    Error: Invalid CG SAM record=260783 67  simulatedSequence1  918 55  23=6N3=1X1D6= simulatedSequence1  1336  418 CGCGGAACTCGTCAGAAATGACGGGTCAGTGTC +]Q[T)/!]!^SX\LN);3UUXP]NV(2TVO3O AS:i:3  MQ:i:255  XU:Z:5=2B1R19=6N3=1X1D6=  XQ:Z:%' XR:Z:C  XA:i:4  IH:i:1  NH:i:1

    final SAMRecord samrec = new SAMRecord(null);
    samrec.setAlignmentStart(3);
    samrec.setCigarString("23=6N3=1X1D6=");
    samrec.setReadString("  CGCGGAACTCGTCAGAAATGACGGGTCAGTGTC".replaceAll(" ", ""));
    final String template = "GTCGCGGAACTCGTCAGAAATGACGACTCTAGGTGGAGTGTCGCAGTCGAGG".replaceAll(" ", "");
    // orig read             ATAAG
    //                        GA
    //                         GT  GAAGGAGTGGCACTTC
    samrec.setBaseQualityString("+]Q[T)/!]!^SX\\LN);3UUXP]NV(2TVO3O");
    samrec.setAttribute(SamUtils.ATTRIBUTE_ALIGNMENT_SCORE, 3);
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "5=2B1R19=6N3=1X1D6=");
    samrec.setAttribute(SamUtils.CG_OVERLAP_QUALITY, "%'");
    samrec.setAttribute(SamUtils.CG_READ_DELTA, "C");
    samrec.setFlags(67);

    try (MemoryPrintStream mps = new MemoryPrintStream()) {
      final SamValidator sv = new SamValidator(mps.printStream(), mps.printStream(), true, false, false, false, false, getParams(false), false);
      final byte[] templ = stringToDna(template);
      final PileUp p = new PileUp(template.length());

      assertEquals(mps.toString(), 3, sv.isAtExpectedRef(templ, samrec, p));
    }
  }

  public void testCgGotohProb3() {
    final SAMRecord samrec = new SAMRecord(null);
    samrec.setAlignmentStart(5);
    samrec.setCigarString("1=1X4=1X17=6N10=");
    samrec.setReadString("       TNAGACNGGTAAAATATGAAGTGAAAGGGAGCTT".replaceAll(" ", ""));
    final String template = "GGATTGAGACTGGTAAAATATGAAGTGACCACCAAAGGGAGCTTGAGAGA".replaceAll(" ", "");
    samrec.setBaseQualityString("+]Q[T)/!]!^SX\\LN);3UUXP]NV(2TVO3O");
    samrec.setAttribute(SamUtils.ATTRIBUTE_ALIGNMENT_SCORE, 1);
    samrec.setAttribute(SamUtils.CG_SUPER_CIGAR, "1=1R3=1B1X1=1R17=6N10=");
    samrec.setAttribute(SamUtils.CG_OVERLAP_QUALITY, "+");
    samrec.setAttribute(SamUtils.CG_READ_DELTA, "T");
    samrec.setFlags(179);

    try (MemoryPrintStream mps = new MemoryPrintStream()) {
      SamValidator sv = new SamValidator(mps.printStream(), mps.printStream(), true, false, false, false, false, getParams(false), false);
      byte[] templ = stringToDna(template);
      PileUp p = new PileUp(template.length());

      assertEquals(mps.toString(), 0, sv.isAtExpectedRef(templ, samrec, p));

      sv = new SamValidator(mps.printStream(), mps.printStream(), true, false, false, false, false, getParams(true), false);
      templ = stringToDna(template);
      p = new PileUp(template.length());

      assertEquals(mps.toString(), 2, sv.isAtExpectedRef(templ, samrec, p));
    }
  }
}
