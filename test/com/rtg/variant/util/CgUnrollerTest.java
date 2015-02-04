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

import static com.rtg.util.StringUtils.LS;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.Locale;

import com.rtg.launcher.GlobalFlags;
import com.rtg.mode.DnaUtils;
import com.rtg.ngs.CgMapCli;
import com.rtg.reader.FastaUtils;
import com.rtg.reader.PrereadArm;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.sam.SamFileAndRecord;
import com.rtg.sam.SamUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.test.FileHelper;
import com.rtg.variant.VariantAlignmentRecord;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;

import junit.framework.TestCase;

/**
 */
public class CgUnrollerTest extends TestCase {

  static final String REF = "TGTTCCAGCACCATTTGTTGAAAAAACTGTCTTTTCTACTGACACTTGTT";

  @Override
  public void setUp() {
    GlobalFlags.resetAccessedStatus();
  }

  public void testBug1296() {
    final byte[] template = DnaUtils.encodeString(REF);
    final SAMRecord rec = new SAMRecord(new SAMFileHeader());
    rec.setReadString("GTTCCAGCACTTGAAAAAACTGTCTTTTC");
    rec.setBaseQualityString("9899441797:<<;;=<<;:::::8:::9");
    rec.setFlags(131);
    rec.setAlignmentStart(2);
    rec.setCigarString("10=6N19=");
    rec.setAttribute(SamUtils.CG_SUPER_CIGAR, "10=6N16=2I2=4B5=");
    rec.setAttribute(SamUtils.CG_OVERLAP_QUALITY, ";54577");
    rec.setAttribute(SamUtils.CG_READ_DELTA, "TC");
    final CgUnroller.OrientedRead read = CgUnroller.unrollCgRead(new VariantAlignmentRecord(rec), template);
    final String readx = new StringBuilder("GTTCCAGCACTTGAAAAAACTGTCTTTCTTTTTTC").reverse().toString();
    final String qual = new StringBuilder("9899441797:<<;;=<<;:::::;545778:::9").reverse().toString();
    final byte[] qual2 = FastaUtils.asciiToRawQuality(qual);

    assertTrue("expected : " + Arrays.toString(FastaUtils.rawToAsciiQuality(qual.getBytes())) + "\nactual : " + Arrays.toString(read.getQuality()),
      Arrays.equals(qual2, read.getQuality()));
    assertTrue("expected : " + Arrays.toString(readx.getBytes()) + "\nactual : " + Arrays.toString(read.getRead()), Arrays.equals(readx.getBytes(), read.getRead()));
  }

  public void testCgUnrollNoOverlap() {
    final SAMRecord rec = new SAMRecord(new SAMFileHeader());
    rec.setReadString("TATTAGGATTGAGACTGGTAAAATGGNCCACCAAG");
    rec.setBaseQualityString("/725361840-525251.68,0,.52!222254/2");
    assertNull(CgUnroller.unrollCgRead(new VariantAlignmentRecord(rec), null));
    rec.setReadPairedFlag(true);
    final CgUnroller.OrientedRead a = CgUnroller.unrollCgRead(new VariantAlignmentRecord(rec), null);
    assertFalse(a.isCgOverlapOnLeft());
    assertEquals("GAACCACCNGGTAAAATGGTCAGAGTTAGGATTAT", new String(a.getRead()));
    assertEquals("2/452222!25.,0,86.152525-048163527/", new String(FastaUtils.rawToAsciiQuality(a.getQuality())));
    rec.setFirstOfPairFlag(true);
    final CgUnroller.OrientedRead b = CgUnroller.unrollCgRead(new VariantAlignmentRecord(rec), null);
    assertTrue(b.isCgOverlapOnLeft());
    assertEquals("TATTAGGATTGAGACTGGTAAAATGGNCCACCAAG", new String(b.getRead()));
    assertEquals("/725361840-525251.68,0,.52!222254/2", new String(FastaUtils.rawToAsciiQuality(b.getQuality())));
    rec.setReadNegativeStrandFlag(true);
    final CgUnroller.OrientedRead c = CgUnroller.unrollCgRead(new VariantAlignmentRecord(rec), null);
    assertFalse(c.isCgOverlapOnLeft());
    assertEquals("GAACCACCNGGTAAAATGGTCAGAGTTAGGATTAT", new String(c.getRead()));
    assertEquals("2/452222!25.,0,86.152525-048163527/", new String(FastaUtils.rawToAsciiQuality(c.getQuality())));
    rec.setFirstOfPairFlag(false);
    final CgUnroller.OrientedRead d = CgUnroller.unrollCgRead(new VariantAlignmentRecord(rec), null);
    assertTrue(d.isCgOverlapOnLeft());
    assertEquals("TATTAGGATTGAGACTGGTAAAATGGNCCACCAAG", new String(d.getRead()));
    assertEquals("/725361840-525251.68,0,.52!222254/2", new String(FastaUtils.rawToAsciiQuality(d.getQuality())));
  }

  public void testCgUnrollNoOverlapWithSuperCigarForNoGoodReason() {
    final byte[] template = DnaUtils.encodeString("TATTAGGATTGAGACTGGTAAAATGNNNNNGNCCACCAAG".replaceAll(" ", ""));
    final SAMRecord rec = new SAMRecord(new SAMFileHeader());
    rec.setReadString("TATTAGGATTGAGACTGGTAAAATGGNCCACCAAG");
    rec.setBaseQualityString("/725361840-525251.68,0,.52!222254/2");
    assertNull(CgUnroller.unrollCgRead(new VariantAlignmentRecord(rec), null));
    rec.setReadPairedFlag(true);
    rec.setAlignmentStart(1);
    rec.setAttribute(SamUtils.CG_SUPER_CIGAR, "25=5N10=");
    final CgUnroller.OrientedRead a = CgUnroller.unrollCgRead(new VariantAlignmentRecord(rec), template);
    assertNotNull(a);      //not convinced that this can't return the expected result, but this case shouldn't happen in the real world anyway.
    assertFalse(a.isCgOverlapOnLeft());
    assertEquals("GAACCACCNGGTAAAATGGTCAGAGTTAGGATTAT", new String(a.getRead()));
    assertEquals("2/452222!25.,0,86.152525-048163527/", new String(FastaUtils.rawToAsciiQuality(a.getQuality())));
    rec.setFirstOfPairFlag(true);
    final CgUnroller.OrientedRead b = CgUnroller.unrollCgRead(new VariantAlignmentRecord(rec), template);
    assertNotNull(b);
    assertTrue(b.isCgOverlapOnLeft());
    assertEquals("TATTAGGATTGAGACTGGTAAAATGGNCCACCAAG", new String(b.getRead()));
    assertEquals("/725361840-525251.68,0,.52!222254/2", new String(FastaUtils.rawToAsciiQuality(b.getQuality())));
    rec.setReadNegativeStrandFlag(true);
    final CgUnroller.OrientedRead c = CgUnroller.unrollCgRead(new VariantAlignmentRecord(rec), template);
    assertNotNull(c);
    assertFalse(c.isCgOverlapOnLeft());
    assertEquals("GAACCACCNGGTAAAATGGTCAGAGTTAGGATTAT", new String(c.getRead()));
    assertEquals("2/452222!25.,0,86.152525-048163527/", new String(FastaUtils.rawToAsciiQuality(c.getQuality())));
    rec.setFirstOfPairFlag(false);
    final CgUnroller.OrientedRead d = CgUnroller.unrollCgRead(new VariantAlignmentRecord(rec), template);
    assertNotNull(d);
    assertTrue(d.isCgOverlapOnLeft());
    assertEquals("TATTAGGATTGAGACTGGTAAAATGGNCCACCAAG", new String(d.getRead()));
    assertEquals("/725361840-525251.68,0,.52!222254/2", new String(FastaUtils.rawToAsciiQuality(d.getQuality())));
  }

  public void testCgUnrollTooLong() {
    final SAMRecord rec = new SAMRecord(new SAMFileHeader());
    rec.setReadString("TATTAGGATTGAGACTGGTAAAATGGNCCACCAAGG");
    rec.setBaseQualityString("/725361840-525251.68,0,.52!222254/23");
    rec.setReadPairedFlag(true);
    rec.setFirstOfPairFlag(true);
    assertNull(CgUnroller.unrollCgRead(new VariantAlignmentRecord(rec), null));
  }

  public void testCgUnrollTooShort() {
    final SAMRecord rec = new SAMRecord(new SAMFileHeader());
    rec.setReadString("TATTAGGATTGAGACTGGTAAAATGGNCCACCAG");
    rec.setBaseQualityString("/725361840-525251.68,0,.52!222254/23");
    rec.setReadPairedFlag(true);
    rec.setFirstOfPairFlag(true);
    assertNull(CgUnroller.unrollCgRead(new VariantAlignmentRecord(rec), null));
  }

  public void testCgUnrollOverlapGSGQInvalid() {
    final SAMRecord rec = new SAMRecord(new SAMFileHeader());
    rec.setReadString("AGGATTGAGACTGGTAAAATATGACACCAAAGGG");
    rec.setBaseQualityString("24/5,,*22-4+-*06524/277.1/1101346.");
    rec.setReadPairedFlag(true);
    rec.setAttribute("GS", "TTG");
    rec.setAttribute("GQ", "1");
    rec.setAttribute("GC", "4S1G29S");
    assertNull(CgUnroller.unrollCgRead(new VariantAlignmentRecord(rec), null));
    rec.setAttribute("GS", "TTGG");
    assertNull(CgUnroller.unrollCgRead(new VariantAlignmentRecord(rec), null));
  }

  public void testCgUnrollOverlapGCInvalid() {
    final SAMRecord rec = new SAMRecord(new SAMFileHeader());
    rec.setReadString("AGGATTGAGACTGGTAAAATATGACACCAAAGG");
//    rec.setBaseQualityString("24/5,,*22-4+-*06524/277.1/1101346.");
    rec.setReadPairedFlag(true);
    rec.setAttribute("GS", "TTGG");
    rec.setAttribute("GQ", "11");
    assertNull(CgUnroller.unrollCgRead(new VariantAlignmentRecord(rec), null));
    rec.setAttribute("GC", "35G");
    assertNull(CgUnroller.unrollCgRead(new VariantAlignmentRecord(rec), null));
  }

  public void testCgUnrollOverlap1() {
    final byte[] template = DnaUtils.encodeString("AGGAT TGAGACTGGTAAAATATGA TTTTTT CACCAAAGGG".replaceAll(" ", ""));
    final SAMRecord rec = new SAMRecord(new SAMFileHeader());
    rec.setReadString("AGGATTGAGACTGGTAAAATATGACACCAAAGGG");
    rec.setBaseQualityString("24/5,,*22-4+-*06524/277.1/1101346.");
    rec.setReadPairedFlag(true);
    rec.setAlignmentStart(1);
    //    rec.setAttribute("GS", "TT");
    //    rec.setAttribute("GQ", "1");
    //    rec.setAttribute("GC", "4S1G29S");
    rec.setAttribute(SamUtils.CG_SUPER_CIGAR, "5=1B20=6N10=");
    // read delta is missing, since this is a perfect match.
    rec.setAttribute("XQ", "1");
    final CgUnroller.OrientedRead a = CgUnroller.unrollCgRead(new VariantAlignmentRecord(rec), template);
    assertFalse(a.isCgOverlapOnLeft());
    assertEquals("GGGAAACCACAGTATAAAATGGTCAGAGTTTAGGA", new String(a.getRead()));
    assertEquals(".6431011/1.772/42560*-+4-22*,1,5/42", new String(FastaUtils.rawToAsciiQuality(a.getQuality())));
    rec.setFirstOfPairFlag(true);
    final CgUnroller.OrientedRead b = CgUnroller.unrollCgRead(new VariantAlignmentRecord(rec), template);
    assertTrue(b.isCgOverlapOnLeft());
    assertEquals("AGGATTTGAGACTGGTAAAATATGACACCAAAGGG", new String(b.getRead()));
    assertEquals("24/5,1,*22-4+-*06524/277.1/1101346.", new String(FastaUtils.rawToAsciiQuality(b.getQuality())));
    rec.setReadNegativeStrandFlag(true);
    final CgUnroller.OrientedRead c = CgUnroller.unrollCgRead(new VariantAlignmentRecord(rec), template);
    assertFalse(c.isCgOverlapOnLeft());
    assertEquals("GGGAAACCACAGTATAAAATGGTCAGAGTTTAGGA", new String(c.getRead()));
    assertEquals(".6431011/1.772/42560*-+4-22*,1,5/42", new String(FastaUtils.rawToAsciiQuality(c.getQuality())));
    rec.setFirstOfPairFlag(false);
    final CgUnroller.OrientedRead d = CgUnroller.unrollCgRead(new VariantAlignmentRecord(rec), template);
    assertTrue(d.isCgOverlapOnLeft());
    assertEquals("AGGATTTGAGACTGGTAAAATATGACACCAAAGGG", new String(d.getRead()));
    assertEquals("24/5,1,*22-4+-*06524/277.1/1101346.", new String(FastaUtils.rawToAsciiQuality(d.getQuality())));
  }

  public void testCgUnrollOverlap2() {
    final byte[] template = DnaUtils.encodeString("GNAAA AGTTATGATTTCAAACGT NNNN CCTTTTTAGA".replaceAll(" ", ""));
    final SAMRecord rec = new SAMRecord(new SAMFileHeader());
    rec.setReadString("GNAAAAGTTATGATTTCAAACGTCCTTTTTAGA");
    rec.setBaseQualityString("1!02/550-561252316+.-.55.,-+#5245");
    rec.setReadPairedFlag(true);
    rec.setAlignmentStart(1);
    //    rec.setAttribute("GS", "AAAC");
    //    rec.setAttribute("GQ", "12");
    //    rec.setAttribute("GC", "3S2G28S");
    rec.setAttribute(SamUtils.CG_SUPER_CIGAR, "5=2B1=1X18=4N10=");
    rec.setAttribute("XR", "C");
    rec.setAttribute("XQ", "12");
    final CgUnroller.OrientedRead a = CgUnroller.unrollCgRead(new VariantAlignmentRecord(rec), template);
    assertFalse(a.isCgOverlapOnLeft());
    assertEquals("AGATTTTTCCTGCAAACTTTAGTATTGACAAAANG", new String(a.getRead()));
    assertEquals("5425#+-,.55.-.+613252165-05521/20!1", new String(FastaUtils.rawToAsciiQuality(a.getQuality())));
    rec.setFirstOfPairFlag(true);
    final CgUnroller.OrientedRead b = CgUnroller.unrollCgRead(new VariantAlignmentRecord(rec), template);
    assertTrue(b.isCgOverlapOnLeft());
    assertEquals("GNAAAACAGTTATGATTTCAAACGTCCTTTTTAGA", new String(b.getRead()));
    assertEquals("1!02/12550-561252316+.-.55.,-+#5245", new String(FastaUtils.rawToAsciiQuality(b.getQuality())));
    rec.setReadNegativeStrandFlag(true);
    final CgUnroller.OrientedRead c = CgUnroller.unrollCgRead(new VariantAlignmentRecord(rec), template);
    assertFalse(c.isCgOverlapOnLeft());
    assertEquals("AGATTTTTCCTGCAAACTTTAGTATTGACAAAANG", new String(c.getRead()));
    assertEquals("5425#+-,.55.-.+613252165-05521/20!1", new String(FastaUtils.rawToAsciiQuality(c.getQuality())));
    rec.setFirstOfPairFlag(false);
    final CgUnroller.OrientedRead d = CgUnroller.unrollCgRead(new VariantAlignmentRecord(rec), template);
    assertTrue(d.isCgOverlapOnLeft());
    assertEquals("GNAAAACAGTTATGATTTCAAACGTCCTTTTTAGA", new String(d.getRead()));
    assertEquals("1!02/12550-561252316+.-.55.,-+#5245", new String(FastaUtils.rawToAsciiQuality(d.getQuality())));
  }

  public void testCgUnrollOverlap2Oldschool() {
    final byte[] template = DnaUtils.encodeString("GNAAA AGTTATGATTTCAAACGT NNNN CCTTTTTAGA".replaceAll(" ", ""));
    final SAMRecord rec = new SAMRecord(new SAMFileHeader());
    rec.setReadString("GNAAAAGTTATGATTTCAAACGTCCTTTTTAGA");
    rec.setBaseQualityString("1!02/550-561252316+.-.55.,-+#5245");
    rec.setReadPairedFlag(true);
    rec.setAlignmentStart(1);
    rec.setAttribute("GS", "AAAC");
    rec.setAttribute("GQ", "12");
    rec.setAttribute("GC", "3S2G28S");
    final CgUnroller.OrientedRead a = CgUnroller.unrollCgRead(new VariantAlignmentRecord(rec), template);
    assertFalse(a.isCgOverlapOnLeft());
    assertEquals("AGATTTTTCCTGCAAACTTTAGTATTGACAAAANG", new String(a.getRead()));
    assertEquals("5425#+-,.55.-.+613252165-05521/20!1", new String(FastaUtils.rawToAsciiQuality(a.getQuality())));
    rec.setFirstOfPairFlag(true);
    final CgUnroller.OrientedRead b = CgUnroller.unrollCgRead(new VariantAlignmentRecord(rec), template);
    assertTrue(b.isCgOverlapOnLeft());
    assertEquals("GNAAAACAGTTATGATTTCAAACGTCCTTTTTAGA", new String(b.getRead()));
    assertEquals("1!02/12550-561252316+.-.55.,-+#5245", new String(FastaUtils.rawToAsciiQuality(b.getQuality())));
    rec.setReadNegativeStrandFlag(true);
    final CgUnroller.OrientedRead c = CgUnroller.unrollCgRead(new VariantAlignmentRecord(rec), template);
    assertFalse(c.isCgOverlapOnLeft());
    assertEquals("AGATTTTTCCTGCAAACTTTAGTATTGACAAAANG", new String(c.getRead()));
    assertEquals("5425#+-,.55.-.+613252165-05521/20!1", new String(FastaUtils.rawToAsciiQuality(c.getQuality())));
    rec.setFirstOfPairFlag(false);
    final CgUnroller.OrientedRead d = CgUnroller.unrollCgRead(new VariantAlignmentRecord(rec), template);
    assertTrue(d.isCgOverlapOnLeft());
    assertEquals("GNAAAACAGTTATGATTTCAAACGTCCTTTTTAGA", new String(d.getRead()));
    assertEquals("1!02/12550-561252316+.-.55.,-+#5245", new String(FastaUtils.rawToAsciiQuality(d.getQuality())));
  }

  public void testCgUnrollNoOverlapNoQuality() {
    final SAMRecord rec = new SAMRecord(new SAMFileHeader());
    rec.setReadString("TATTAGGATTGAGACTGGTAAAATGGNCCACCAAG");
    assertNull(CgUnroller.unrollCgRead(new VariantAlignmentRecord(rec), null));
    rec.setReadPairedFlag(true);
    final CgUnroller.OrientedRead a = CgUnroller.unrollCgRead(new VariantAlignmentRecord(rec), null);
    assertFalse(a.isCgOverlapOnLeft());
    assertEquals("GAACCACCNGGTAAAATGGTCAGAGTTAGGATTAT", new String(a.getRead()));
    assertNull(a.getQuality());
    rec.setFirstOfPairFlag(true);
    final CgUnroller.OrientedRead b = CgUnroller.unrollCgRead(new VariantAlignmentRecord(rec), null);
    assertTrue(b.isCgOverlapOnLeft());
    assertEquals("TATTAGGATTGAGACTGGTAAAATGGNCCACCAAG", new String(b.getRead()));
    assertNull(b.getQuality());
  }

  public void testCgUnrollOverlap2NoQuality() {
    final byte[] template = DnaUtils.encodeString("GNAAA AGTTATGATTTCAAACGT NNNNNNNN CCTTTTTAGA".replaceAll(" ", ""));
    final SAMRecord rec = new SAMRecord(new SAMFileHeader());
    rec.setReadString("GNAAAAGTTATGATTTCAAACGTCCTTTTTAGA");
    rec.setAlignmentStart(1);
    //    rec.setAttribute("GS", "AAAC");
    //    rec.setAttribute("GQ", "12");
    //    rec.setAttribute("GC", "3S2G28S");
    rec.setAttribute(SamUtils.CG_SUPER_CIGAR, "5=2B1=1X18=8N10=");
    rec.setAttribute("XR", "C");
    rec.setReadPairedFlag(true);
    rec.setFirstOfPairFlag(true);
    rec.setReadNegativeStrandFlag(true);
    final CgUnroller.OrientedRead c = CgUnroller.unrollCgRead(new VariantAlignmentRecord(rec), template);
    assertFalse(c.isCgOverlapOnLeft());
    assertEquals("AGATTTTTCCTGCAAACTTTAGTATTGACAAAANG", new String(c.getRead()));
    assertNull(c.getQuality());
  }

  public void testCgBadSuperCigar() {
    final ByteArrayOutputStream bos = new ByteArrayOutputStream();
    Diagnostic.setLogStream(new PrintStream(bos));
    final byte[] template = DnaUtils.encodeString("GNAAA AGTTATGATTTCAAACGT NNNNNNNN CCTTTTTAGA".replaceAll(" ", ""));
    final SAMRecord rec = new SAMRecord(new SAMFileHeader());
    //rec.setCgReadString("GNAAAAGTTATGATTTCAAACGTCCTTTTTAGA");
    rec.setAlignmentStart(1);
    rec.setAttribute(SamUtils.CG_SUPER_CIGAR, "5=2B1=1X18=8N10=");
    rec.setAttribute("XR", "");  // too short, should be length 1
    rec.setReadPairedFlag(true);
    rec.setFirstOfPairFlag(true);
    rec.setReadNegativeStrandFlag(true);
    final CgUnroller.OrientedRead c = CgUnroller.unrollCgRead(new VariantAlignmentRecord(rec), template);
    assertNull(c);
    Diagnostic.setLogStream(); // closes bos stream.
    assertTrue(bos.toString().contains("Ignored SAM CG record due to Ill-formed cigar/read delta: 5=2B1=1X18=8N10=/"));
  }

  public void testUnrollCgAlan() {
    final byte[] tmpl = DnaUtils.encodeString("TTCGAGTATGGAGCCAGTCCTCATAAGTCGGCCCCTTCCTATGTCC");
    final SAMFileHeader sft3 = new SAMFileHeader();
    final SAMRecord rec = new SAMRecord(sft3);
    rec.setReadName("blah");
    rec.setReadString("CGAGTATGGACCTCATAAGTCGGCCCCTGGT");       //used just for the length
    rec.setReferenceName("sdjr");
    rec.setFlags(131);
    rec.setAlignmentStart(3);
    rec.setInferredInsertSize(0); // Default, overwritten for valid paired-end
    rec.setMateAlignmentStart(0); // Default, overwritten for valid paired-end
    rec.setMateReferenceName("*");
    rec.setAttribute(SamUtils.CG_SUPER_CIGAR, "10=6N20=4B2=2I1=");
    rec.setAttribute(SamUtils.CG_READ_DELTA, "GG");

    final String expRead = "TGGTCCTTCCCCGGCTGAATACTCCAGGTATGAGC";
    final CgUnroller.OrientedRead or = CgUnroller.unrollCgRead(new VariantAlignmentRecord(rec), tmpl);
    assertEquals(expRead.length(), or.getRead().length);
    for (int i = 0; i < expRead.length(); i++) {
      assertEquals(expRead.charAt(i), or.getRead()[i]);
    }
  }

//Test that the conventions for CG SAM representation are the same in SAM writing and unrolling

  private static final String CG_QUALITY = "+" + LS + "00000" + "0000000000" + "0000000000" + "0000000000" + LS;

  private static final String READ_L_F = ""
    + "acgTA" + "TAgtacgact" + "acgtacgtac" + "acgtacgtac"
    ;
  private static final String READ_R_F = ""
    + "ttccaaggct" + "ttccaaggct" + "gcatatttGC" + "GCtac"
    ;

  public void testEndToEndOverlap2SameF() throws IOException {
    checkCGOverlap(READ_L_F, READ_R_F, 2, false, "23=6N10=/10=6N23=");
  }

  public void testEndToEndOverlap2SameR() throws IOException {
    checkCGOverlap(READ_L_F, READ_R_F, 2, true, "23=6N10=/10=6N23=");
  }

  private static final String READ_L_F_DIFF = ""
    + "acgTA" + "GAgtacgact" + "acgtacgtac" + "acgtacgtac"
    ;
  private static final String READ_R_F_DIFF = ""
    + "ttccaaggct" + "ttccaaggct" + "gcatatttGC" + "ACtac"
    ;

  public void testEndToEndOverlap2DiffF() throws IOException {
    //checkCGOverlap(READ_L_F_DIFF, READ_R_F_DIFF, 2, false, "4=1I19=6N10=/10=6N18=1X4=");
    // with CG gotoh aligner:
    checkCGOverlap(READ_L_F_DIFF, READ_R_F_DIFF, 2, false, "23=6N10=/10=6N23=");
  }

  public void testEndToEndOverlap2DiffR() throws IOException {
    //checkCGOverlap(READ_L_F_DIFF, READ_R_F_DIFF, 2, true, "4=1I19=6N10=/10=6N18=1X4=");
    // with CG gotoh aligner:
    checkCGOverlap(READ_L_F_DIFF, READ_R_F_DIFF, 2, true, "23=6N10=/10=6N23=");
  }

  public void testCGMapbug1137() throws IOException {
    final String leftStr = READ_L_F.substring(0, 5) + READ_L_F.substring(5 + 2, 25) + "ACTGAC" + READ_L_F.substring(25);
    assertEquals(35 + 6 - 2, leftStr.length());
    final String rightStr = READ_R_F.substring(0, 10) + "ACTGAC" + READ_R_F.substring(10, 30 - 2) + READ_R_F.substring(30);
    assertEquals(35 + 6 - 2, rightStr.length());
    final String templateSeq = ""
      + ">left" + LS + leftStr + LS
      + ">right" + LS + rightStr + LS
      ;
    //System.err.println(templateSeq);
    final PrintStream err = System.err;
    final File parent = FileUtils.createTempDir("testSamUtils", "overlap2Same");
    try {
      System.setErr(TestUtils.getNullPrintStream());
      final File reads = new File(parent, "reads");
      final File left = new File(reads, "left");
      final File right = new File(reads, "right");
      final String left0 = READ_L_F;
      ReaderTestUtils.getReaderDNAFastqCG("@reL" + LS + left0 + CG_QUALITY, left, PrereadArm.LEFT);
      final String right0 = READ_R_F;
      ReaderTestUtils.getReaderDNAFastqCG("@reR" + LS + right0 + CG_QUALITY, right, PrereadArm.RIGHT);

      final File template = new File(parent, "template");
      ReaderTestUtils.getReaderDNA(templateSeq, template, null);
      final File output = new File(parent, "out");
      final CgMapCli cgmap = new CgMapCli();
      assertEquals(0, cgmap.mainInit(new String[] {"-o", output.getPath(), "-i", reads.getPath(), "-t", template.getPath(), "-Z", "-T", "1", "--sam", "--no-merge"}, /*err*/TestUtils.getNullOutputStream(), System.err));
      // The following lines were failing due to a bug
      final String str = unrollCGFromSam(new File(output, "unmated.sam"), leftStr, rightStr);
      final String exp = (READ_L_F + LS + READ_R_F + LS).toUpperCase(Locale.getDefault());
      assertEquals(exp, str);

      //System.err.println(str);
    } finally {
      System.setErr(err);
      //System.err.println(parent.getPath());
      FileHelper.deleteAll(parent);
    }
  }

  private void checkCGOverlap(final String leftRead, final String rightRead, final int overlap, final boolean reverse, final String cigars) throws IOException {
    Diagnostic.setLogStream(new PrintStream(new ByteArrayOutputStream()));
    //NOTE: with the old CG aligner, these tests passed regardless of which side of the overlap we discarded!
    //      (because the GS SAM attribute stored both sides of the overlap).
    //But with the CG Gotoh aligner, less is stored in the SAM file, so this end-to-end test works
    //      only if the middle-most side of the overlap region is discarded.
    final String leftStr = leftRead.substring(0, 5) + leftRead.substring(5 + overlap, 25) + "ACTGAC" + leftRead.substring(25);
    //final String leftStr = leftRead.substring(0, 5 - overlap) + leftRead.substring(5, 25) + "ACTGAC" + leftRead.substring(25);
    assertEquals(35 + 6 - overlap, leftStr.length());
    final String rightStr = rightRead.substring(0, 10) + "ACTGAC" + rightRead.substring(10, 30 - overlap) + rightRead.substring(30);
    //final String rightStr = rightRead.substring(0, 10) + "ACTGAC" + rightRead.substring(10, 30) + rightRead.substring(30 + overlap);
    assertEquals(35 + 6 - overlap, rightStr.length());
    final String templateSeq = ""
      + ">left" + LS + leftStr + LS
      + ">right" + LS + rightStr + LS
      ;
    //    System.out.println("Template:  " + templateSeq.replaceAll(LS, " "));
    final MemoryPrintStream err = new MemoryPrintStream();
    final File parent = FileUtils.createTempDir("testSamUtils", "overlap2Same");
    try {
      final File reads = new File(parent, "reads");
      final File left = new File(reads, "left");
      final File right = new File(reads, "right");
      final String left0 = reverse ? DnaUtils.reverseComplement(rightRead) : leftRead;
      ReaderTestUtils.getReaderDNAFastqCG("@reL" + LS + left0 + CG_QUALITY, left, PrereadArm.LEFT);
      final String right0 = reverse ? DnaUtils.reverseComplement(leftRead) : rightRead;
      ReaderTestUtils.getReaderDNAFastqCG("@reR" + LS + right0 + CG_QUALITY, right, PrereadArm.RIGHT);

      final File template = new File(parent, "template");
      ReaderTestUtils.getReaderDNA(templateSeq, template, null);
      final File output = new File(parent, "out");
      final CgMapCli cgmap = new CgMapCli();
      final int code = cgmap.mainInit(new String[] {"-o", output.getPath(), "-i", reads.getPath(), "-t", template.getPath(), "-Z", "-T", "1", "--sam", "--no-merge"}, err.outputStream(), err.printStream());
      assertEquals(err.toString(), 0, code);
      final String str = unrollCGFromSam(new File(output, "unmated.sam"), leftStr, rightStr);
      final String exp = (leftRead + LS + rightRead + LS).toUpperCase(Locale.getDefault());
      //      System.out.println("Expecting:     " + exp.replaceAll(LS, "                     "));
      //      System.out.println("UnrollFromSam: " + str.replaceAll(LS, "                     "));
      assertEquals(exp, str);
      assertEquals(cigars, mCigarStrings);
    } finally {
      //System.err.println(parent.getPath());
      FileHelper.deleteAll(parent);
    }
  }

  private String mCigarStrings; // stores the cigar found by the most recent call to unrollCGFromSam.

  private String unrollCGFromSam(final File samFile, String leftTemplate, String rightTemplate) throws IOException {
    final StringBuilder sb = new StringBuilder();
    final StringBuilder sbCigar = new StringBuilder();
    final SamFileAndRecord iter = new SamFileAndRecord(samFile, 0);
    while (iter.hasNext()) {
      final SAMRecord sam = iter.next();
      final byte[] template = DnaUtils.encodeString(sam.getReferenceName().equals("left") ? leftTemplate : rightTemplate);
      sbCigar.append(sam.getCigarString()).append("/");
      final VariantAlignmentRecord rec = new VariantAlignmentRecord(sam);
      final CgUnroller.OrientedRead unr = CgUnroller.unrollCgRead(rec, template);
      assertNotNull(unr);
      final byte[] ur = unr.getRead();
      final String ustr = new String(ur);
      if (unr.isCgOverlapOnLeft()) {
        sb.append(ustr).append(LS);
      } else {
        //right arm reverse it
        final StringBuilder sbr = new StringBuilder();
        sbr.append(ustr);
        sb.append(sbr.reverse()).append(LS);
      }
    }
    iter.close();
    sbCigar.deleteCharAt(sbCigar.length() - 1);
    mCigarStrings = sbCigar.toString();
    return sb.toString();
  }

  public void testNoOverlapUnroll() {
    //Error: Invalid CG SAM record=586345 115 simulatedSequence5  27460 55  10=6N4=1X2=1X17=  simulatedSequence5  27059 -401  CAGAANGGAANCTCGGTCTCAGTTGTCCTCAGCAA ;+I.!-_[O<18IA^J/7+/<*EX`OY2X?YG$BB AS:i:2  MQ:i:255  XU:Z:5=1R4=6N1R3=1X2=1X17=  XR:Z:GC XA:i:3  IH:i:1  NH:i:1

    final byte[] tmpl = DnaUtils.encodeString("CAGAANGGAATTTTTTNCTCAGTATCAGTTGTCCTCAGCAA");
    final SAMFileHeader sft3 = new SAMFileHeader();
    final SAMRecord rec = new SAMRecord(sft3);
    rec.setReadName("blah");
    rec.setReadString("CAGAANGGAANCTCGGTCTCAGTTGTCCTCAGCAA");       //used just for the length
    rec.setReferenceName("sdjr");
    rec.setFlags(115);
    rec.setAlignmentStart(1);
    rec.setInferredInsertSize(0); // Default, overwritten for valid paired-end
    rec.setMateAlignmentStart(0); // Default, overwritten for valid paired-end
    rec.setMateReferenceName("*");
    rec.setAttribute(SamUtils.CG_SUPER_CIGAR, "5=1R4=6N1R3=1X2=1X17=");
    rec.setAttribute(SamUtils.CG_READ_DELTA, "GC");
    rec.setBaseQualityString(";+I.!-_[O<18IA^J/7+/<*EX`OY2X?YG$BB");

    final String expRead = "AACGACTCCTGTTGACTCTGGCTCNAAGGNAAGAC";
    final String expQual = "BB$GY?X2YO`XE*</+7/J^AI81<O[_-!.I+;";
    final CgUnroller.OrientedRead or = CgUnroller.unrollCgRead(new VariantAlignmentRecord(rec), tmpl);
    assertNotNull(or);
    assertEquals(expRead.length(), or.getRead().length);
    for (int i = 0; i < expRead.length(); i++) {
      assertEquals(expRead.charAt(i), or.getRead()[i]);
    }
    assertEquals(expQual.length(), or.getQuality().length);
    for (int i = 0; i < expQual.length(); i++) {
      assertEquals(expQual.charAt(i), FastaUtils.rawToAsciiQuality(or.getQuality()[i]));
    }
  }

}
