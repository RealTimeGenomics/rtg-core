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
package com.rtg.alignment;

import static com.rtg.util.StringUtils.LS;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Locale;

import com.rtg.launcher.AbstractNanoTest;
import com.rtg.mode.DnaUtils;
import com.rtg.ngs.tempstage.BinaryTempFileRecord;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.PortableRandom;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.variant.MachineErrorParams;
import com.rtg.variant.MachineErrorParamsBuilder;
import com.rtg.variant.realign.RealignParams;
import com.rtg.variant.realign.RealignParamsImplementation;
import com.rtg.variant.realign.AbstractScoreMatrixCGTest;

/**
 */
public class CgGotohEditDistanceTest extends AbstractNanoTest {

  private RealignParams mParams;
  private CgGotohEditDistance mCgGotoh;
  private String mLog;
  private AlignmentResult mAlignment;
  private BinaryTempFileRecord mSamRecord;

  @Override
  public void tearDown() throws IOException {
    mParams = null;
    mCgGotoh = null;
    mLog = null;
    mAlignment = null;
    mSamRecord = null;
    super.tearDown();
  }

  protected BidirectionalEditDistance getEditDistanceInstance(int unknownsPenalty, String errorsFile, boolean v2) {
    try {
      mParams = new RealignParamsImplementation(new MachineErrorParamsBuilder().errors(errorsFile).create());
    } catch (final InvalidParamsException | IOException e) {
      throw new RuntimeException(e);
    }
    mCgGotoh = new CgGotohEditDistance(7, mParams, unknownsPenalty, v2);
    return new RcEditDistance(mCgGotoh);
  }

  protected BidirectionalEditDistance getEditDistanceInstance(int unknownsPenalty) {
    return getEditDistanceInstance(unknownsPenalty, "cg_test_errors-080412", false);
  }

  protected void checkCG(String read, String template, String x1, String x2, String m, String xCigar, int score, boolean rc, boolean left) {
    checkCG(read, template, x1, x2, m, xCigar, score, true, rc, 0, left);
  }

  protected void checkCG(String read, String template, String x1, String x2,
      String m, String xCigar, int score, boolean first, boolean rc, int start, boolean left) {
    checkCG(read, template, x1, x2, m, xCigar, score, first, rc, start, left, 0, false);
  }

  protected void checkCG(String read, String template, String x1, String x2,
      String actions, String xCigar, int score, boolean first, boolean rc, int start, boolean left, int unknownsPenalty, boolean v2) {
    final byte[] s1 = DnaUtils.encodeString(read.replaceAll(" ", "").toLowerCase(Locale.ROOT));
    final byte[] s2 = DnaUtils.encodeString(template.replaceAll(" ", "").toLowerCase(Locale.ROOT));
    final ByteArrayOutputStream bos = new ByteArrayOutputStream();
    Diagnostic.setLogStream(new PrintStream(bos));
    final BidirectionalEditDistance f = getEditDistanceInstance(unknownsPenalty, "cg_test_errors-080412", v2);
    final int[] v = f.calculateEditDistance(s1, s1.length, s2, start, rc, Integer.MAX_VALUE, 7, left);
    f.logStats();
        //System.out.println(mCgGotoh.toString());
        //System.out.println("start: " + ActionsHelper.zeroBasedTemplateStart(v) + " diff:" + (start - ActionsHelper.zeroBasedTemplateStart(v))
         //   + " score: " + ActionsHelper.alignmentScore(v));
        //System.out.println("actions: " + ActionsHelper.toString(v));
    mAlignment = new AlignmentResult(s1, v, s2);
    mAlignment.setIdentifyingInfo(first, rc);

    final String cigar = mAlignment.getCigarString(rc, false);
    //     System.out.println("Ex:" + xCigar + "|" + cigar);
    final String ss = mAlignment.toString();
    //    System.err.println("alignment result: " + ss);
    mAlignment.setRemainingOutput(1, 0);
    mSamRecord = mAlignment.toRecord(false, null, 0, false, false);
    assertEquals(xCigar, new String(mSamRecord.getCigarString()));
    assertEquals(ss, score, ActionsHelper.alignmentScore(v));
    assertEquals(xCigar, cigar);

    assertEquals(x1.replaceAll("\\.", ""), mAlignment.readString());
    assertEquals(actions, ActionsHelper.toString(v));

//    assertEquals(ss, x1 + "\t" + x2 + "\t" + actions, mAlignment.tabularString());
    Diagnostic.setLogStream();
    mLog = bos.toString();
    //System.err.println("log: " + mLog);

    if (!v2) {
      // Since v1 reads have symmetric layout, we cancheck that we get roughly the same alignment if we flip read+template onto the opposite arm.
      // (this is not always true, but the scores should always be pretty close)
      final byte[] s1R = DnaUtils.encodeString(StringUtils.reverse(read).replaceAll(" ", "").toLowerCase(Locale.ROOT));
      final byte[] s2R = DnaUtils.encodeString(StringUtils.reverse(template).replaceAll(" ", "").toLowerCase(Locale.ROOT));
      assertEquals(s1.length, s1R.length);
      assertEquals(s2.length, s2R.length);
      final BidirectionalEditDistance fR = getEditDistanceInstance(unknownsPenalty, "cg_test_errors-080412", v2); // because f.logStats() kills it!
      final int startR = s2.length - (start + s1.length) - CgGotohEditDistance.CG_INSERT_REGION_DEFAULT_SIZE;
      final int[] vR = fR.calculateEditDistance(s1R, s1R.length, s2R, startR, rc, Integer.MAX_VALUE, 7, !left);
//      System.out.println(mCgGotoh.toString());
//      System.out.println("actions: " + ActionsHelper.toString(vR));
//      System.err.println("startR: " + startR);
      assertEquals(ss, score, ActionsHelper.alignmentScore(vR));
      // the actions strings should be nearly the same, but can have a short reversed segment due to an ambiguous insert/delete position.
      final String actsLeft = ActionsHelper.toString(v);
      final String actsLeftRev = StringUtils.reverse(actsLeft);
      final String acts = ActionsHelper.toString(vR);
      assertEquals(actsLeftRev.length(), acts.length());
      assertEquals(actsLeftRev, acts);
    }
  }

  // this is the example in com/rtg/variant/realign/scorematrixtestCG.xls
  public void testCGSpreadsheet() throws IOException {
    final CgGotohEditDistance ed = new CgGotohEditDistance(7, new AbstractScoreMatrixCGTest.MockRealignParamsCG(), 0);
    final String read = "        ATAAA AAGGCGACAT GCCAATGTGT        TTCAACTTTC";
    final String tmpl = "gggggggataAAA   AAGGCGACAT GCCAATGTGT CGCCTTTTTCAACTTTCCGATTAA";
    final byte[] s1 = DnaUtils.encodeString(read.replaceAll(" ", "").toLowerCase(Locale.ROOT));
    final byte[] s2 = DnaUtils.encodeString(tmpl.replaceAll(" ", "").toLowerCase(Locale.ROOT));
    //    final long start = System.nanoTime();
    final int[] v = ed.calculateEditDistance(s1, s1.length, s2, 7, Integer.MAX_VALUE, 7, true);
    //    final long end = System.nanoTime();
    //    System.out.println("time taken in seconds = " + (end - start) / 1.0E9);
    //    System.out.println(ed.toString());
    //    System.out.println("start: " + ActionsHelper.zeroBasedTemplateStart(v) + " score: " + ActionsHelper.alignmentScore(v));
    //    System.out.println("actions: " + ActionsHelper.toString(v));
    final AlignmentResult alignment = new AlignmentResult(s1, v, s2);
    assertEquals("1=1X21=7N10=", alignment.getCigarString(false, false));
    assertEquals(1, ActionsHelper.alignmentScore(v));
    assertEquals(10, ActionsHelper.zeroBasedTemplateStart(v));
    assertEquals(10, alignment.getStart());
    mNano.check("cg-gotoh-edit-distance-spreadsheet-left.txt", ed.toString());
    mNano.check("cg-gotoh-edit-distance-spreadsheet-strip.txt", stripTopBottom(ed.toString())); // Note, same as next test below
  }

  public void testCGSpreadsheetRightArm() throws IOException {
    final CgGotohEditDistance ed = new CgGotohEditDistance(7, new AbstractScoreMatrixCGTest.MockRealignParamsCG(), 0);
    //final String read = "          ATAAA AAGGCGACAT GCCAATGTGT .....  TTCAACTTTC";
    //final String tmpl = "gggggggataAAA   AAGGCGACAT GCCAATGTGT CGCCTTTTTCAACTTTCCGATTAA";
    final String read = "       CTTTCAACTT         TGTGTAACCG TACAGCGGAA AAATA";
    final String tmpl = "AATTAGCCTTTCAACTT TTTCCGC TGTGTAACCG TACAGCGGAA   AAAataggggggg";
    final byte[] s1 = DnaUtils.encodeString(read.replaceAll(" ", "").toLowerCase(Locale.ROOT));
    final byte[] s2 = DnaUtils.encodeString(tmpl.replaceAll(" ", "").toLowerCase(Locale.ROOT));
    final int[] v = ed.calculateEditDistance(s1, s1.length, s2, 11, Integer.MAX_VALUE, 7, false);
    //    System.out.println(ed.toString());
    //    System.out.println("start: " + ActionsHelper.zeroBasedTemplateStart(v) + " score: " + ActionsHelper.alignmentScore(v));
    //    System.out.println("actions: " + ActionsHelper.toString(v));
    final AlignmentResult alignment = new AlignmentResult(s1, v, s2);
    assertEquals("10=7N21=1X1=", alignment.getCigarString(false, false));
    assertEquals(1, ActionsHelper.alignmentScore(v));
    assertEquals(7, ActionsHelper.zeroBasedTemplateStart(v));
    assertEquals(7, alignment.getStart());
    mNano.check("cg-gotoh-edit-distance-spreadsheet-strip.txt", stripTopBottom(ed.toString())); // Note, same as test above
    mNano.check("cg-gotoh-edit-distance-spreadsheet-right.txt", ed.toString());
  }

  // Strip the template numbering in the top and bottom lines.
  public String stripTopBottom(String edStr) {
    final String[] actualLines = edStr.split(LS);
    final StringBuilder sb = new StringBuilder();
    for (int i = 1; i < actualLines.length - 1; ++i) {
      sb.append(actualLines[i]).append(LS);
    }
    return sb.toString();
  }

  private void logContains(final String gap, final String range, final int... values) {
    final StringBuilder sb = new StringBuilder();
    sb.append(gap);
    sb.append(":\t");
    sb.append(range);
    for (final int val : values) {
      sb.append("\t");
      sb.append(val);
    }
    final String expected = sb.toString();
    assertTrue("Missing from log: " + expected + "\nLog: " + mLog, mLog.contains(expected));
  }

  public void testCGSnpInOverlap1() {
    checkCG("        ATAAG ACGGCGACAT GCCAATGTGT        TTCAACTTTC",
            "gggggggataAAA   ACGGCGACAT GCCAATGTGT CGCCTTTTTCAACTTTCCGATTAA",
            "ATAAGGGCGACATGCCAATGTGT.......TTCAACTTTC",
            "AAAACGGCGACATGCCAATGTGTCGCCTTTTTCAACTTTC",
            "=X==XBB====================NNNNNNN==========",
            "1=1X2=1X18=7N10=", 2, true, false, 7, true);
    //System.out.println(mLog);
    assertTrue(mLog.contains("CgGotohEditDistance mismatch=0.036180 delopen=0.001000 delextend=0.181694 insopen=0.002000 insextend=0.122598"));
    logContains("LeftOverlap", "-4..0", 0, 0, 1, 0, 0);
    logContains("LeftSmall", "0..3", 1, 0, 0, 0);
    logContains("LeftLarge", "4..8", 0, 0, 0, 1, 0);
    logContains("RightOverlap", "-4..0", 0, 0, 0, 0, 0);
    logContains("RightSmall", "0..3", 0, 0, 0, 0);
    logContains("RightLarge", "4..8", 0, 0, 0, 0, 0);
  }

  public void testCGSingleSnp() {
    checkCG("        AAAAC ACGGCGACAT GCCAtTGTGT        TTCAACTTTC",
        "gggggggattAAAAC   GGCGACAT GCCAATGTGT CGCCTT TTCAACTTTCCGATTAA",
        "AAAACGGCGACATGCCATTGTGT......TTCAACTTTC",
        "AAAACGGCGACATGCCAATGTGTCGCCTTTTCAACTTTC",
        "=====BB==============X=====NNNNNN==========",
        "17=1X5=6N10=", 1, true, false, 7, true);
    assertEquals(11, mSamRecord.getStartPosition());
  }
  public void testCGSingleSnpShift() {
    checkCG("        AAAAC ACGGCGACAT GCCAtTGTGT        TTCAACTTTC",
        "gggggggattAAAAC   GGCGACAT GCCAATGTGT CGCCTT TTCAACTTTCCGATTAA",
        "AAAACGGCGACATGCCATTGTGT......TTCAACTTTC",
        "AAAACGGCGACATGCCAATGTGTCGCCTTTTCAACTTTC",
        "=====BB==============X=====NNNNNN==========",
        "17=1X5=6N10=", 1, true, false, 5, true);
    assertEquals(11, mSamRecord.getStartPosition());
  }

  public void testCGLowestScore() {
    checkCG("        AAAAC ACGGCGACAT GCCAATGTGT        TTCAACTTTC",
        "gggggggattAAAAC   GGCGACAT GCCAATGTGT CGCCTT TTCAACTTTCCGATTAA",
        "AAAACGGCGACATGCCAATGTGT......TTCAACTTTC",
        "AAAACGGCGACATGCCAATGTGTCGCCTTTTCAACTTTC",
        "=====BB====================NNNNNN==========",
        "23=6N10=", 0, true, false, 7, true);
    assertEquals(11, mSamRecord.getStartPosition());
  }

  public void testCGSnpInOverlap2() {
    checkCG("        ATAAC AGGGCGACAT GCCAATGTGT        TTCAACTTTC",
          "gggggggataAAAAC   GGCGACAT GCCAATGTGT CGCCTTTTTCAACTTTCCGATTAA",
          "ATAACGGCGACATGCCAATGTGT.......TTCAACTTTC",
          "AAAACGGCGACATGCCAATGTGTCGCCTTTTTCAACTTTC",
          "=X===BB=X==================NNNNNNN==========",
          "1=1X21=7N10=", 2, true, false, 11, true);
    assertEquals(10, mAlignment.getStart());
  }

  public void testCGDelete1Mid() {
    checkCG("        ATAAA AAGGCGACAT GCCA TGTGTC        TTCAACTTTC",
            "gggggggataAAA   AAGGCGACAT GCCAATGTGTC GCCTTT TTCAACTTTCCGATTAA",
            "ATAAAGGCGACATGCCA TGTGTC......TTCAACTTTC".replaceAll(" ", ""),
            "AAAAAGGCGACATGCCAATGTGTCGCCTTTTTCAACTTTC",
            "=X===BB==============D======NNNNNN==========",
//            "| ||||||||||||||| ||||||      ||||||||||",
            "1=1X15=1D6=6N10=", 3, true, false, 7, true);
    assertEquals(10, mAlignment.getStart());
  }

  public void testCGDelete4Mid() {
    checkCG("        ATAAA AAGGCGACAT GCCA    GTCGCC        TTCAACTTTC",
            "gggggggataAAA   AAGGCGACAT GCCAATGTGTCGCC TTTTTT TTCAACTTTCCGATTAA",
            "ATAAAGGCGACATGCCA    GTCGCC......TTCAACTTTC".replaceAll(" ", ""),
            "AAAAAGGCGACATGCCAATGTGTCGCCTTTTTTTTCAACTTTC",
            "=X===BB==============DDDD======NNNNNN==========",
            "1=1X15=4D6=6N10=", 6, true, false, 7, true);
  }

  public void testCGInsert1Mid() {
    checkCG("        ATAAA AAGGCGACAT CCAATGTGTC        TTCAACTTTC",
        "GGGGGGGATAAAA   AAGGCGACAT CCA TGTGTC GCCTTT TTCAACTTTCCGATTAA",
        "ATAAAGGCGACATCCAATGTGTC......TTCAACTTTC",
        "AAAAAGGCGACATCCA-TGTGTCGCCTTTTTCAACTTTC",
        "=X===BB=============I======NNNNNN==========",
        "1=1X14=1I6=6N10=", 3, true, false, 7, true);
  }

  public void testCGInsert4Mid() {
    checkCG("        ATAAA AAGGCGACAT GATGTCTCGC        TTCAACTTTC",
      "gggggggataAAA   AAGGCGACAT     TCTCGC TTTTTT TTCAACTTTCCGATTAA",
      "ATAAAGGCGACATGATGTCTCGC......TTCAACTTTC",
      "AAAAAGGCGACAT----TCTCGCTTTTTTTTCAACTTTC",
      "=X===BB==========IIII======NNNNNN==========",
      "1=1X11=4I6=6N10=", 6, true, false, 7, true);
  }

  // overlap of 0 (cg_test_errors priors do allow this!)
  public void testCGOverlap0() {
    checkCG("acgtcacgtcacgtcacgtcacgtc     acgtcacgtc",
            "acgtcacgtcacgtcacgtcacgtcgggggacgtcacgtc",
            "ACGTCACGTCACGTCACGTCACGTC.....ACGTCACGTC",
            "ACGTCACGTCACGTCACGTCACGTCGGGGGACGTCACGTC",
            "=========================NNNNN==========",
            "25=5N10=", 0, false, true);
    assertEquals(0, mAlignment.getStart());
  }
  public void testCGOverlap0Rev() {
    checkCG("acgtcacgtcacgtcacgtcacgtc     acgtcacgtc",
      "gacgtgacgtcccccgacgtgacgtgacgtgacgtgacgt",
      "ACGTCACGTCACGTCACGTCACGTC.....ACGTCACGTC",
      "ACGTCACGTCACGTCACGTCACGTCGGGGGACGTCACGTC",
      "=========================NNNNN==========",
      "10=5N25=", 0, true, true);
    assertEquals(0, mAlignment.getStart());
  }

  // overlap of 4 (cg_test_errors priors do allow this!)
  public void testCGOverlap4() {
    checkCG("acgtc cgtc aacgtcacgtcacgtc     acgtcacgtc",
      "acgtc      aacgtcacgtcacgtcgggggacgtcacgtc",
      "ACGTCAACGTCACGTCACGTC.....ACGTCACGTC",
      "ACGTCAACGTCACGTCACGTCGGGGGACGTCACGTC",
      "=====BBBB====================NNNNN==========",
      "21=5N10=", 0, false, true);
    assertEquals(1, mSamRecord.getStartPosition());
    assertEquals(0, mAlignment.getStart());
  }

  // large gap of 8, small gap of 3, and overlap of 4 (right arm)
  public void testCG834RightArm() {
    checkCG("aacgtcacgt        cacgtcacgt   cacgtc cgtc cgtca",
            "aacgtcacgtggggggggcacgtcacgtgggcacgtc      cgtca",
            "AACGTCACGT........CACGTCACGT...CACGTCCGTCA",   //TODO: fix this: first 10 are "ccccccccccc"
            "AACGTCACGTGGGGGGGGCACGTCACGTGGGCACGTCCGTCA",
            "==========NNNNNNNN==========NNN==========BBBB=====",
            "10=8N10=3N11=", 0, false, false, 0, false);
    logContains("LeftOverlap", "-4..0", 0, 0, 0, 0, 0);
    logContains("LeftSmall", "0..3", 0, 0, 0, 0);
    logContains("LeftLarge", "4..8", 0, 0, 0, 0, 0);
    logContains("RightOverlap", "-4..0", 1, 0, 0, 0, 0);
    logContains("RightSmall", "0..3", 0, 0, 0, 1);
    logContains("RightLarge", "4..8", 0, 0, 0, 0, 1);
    assertEquals(1, mSamRecord.getStartPosition());
  }

  // large gap of 4, small gap of 0, and overlap of 0 (right arm)
  public void testCG400RightArm() {
    checkCG("aacgtcacgt     cacgtcacgtcacgtc cgtc cgtca",
            "aacgtcacgtgggg cacgtcacgtcacgtc cgtc cgtca",
            "AACGTCACGT....CACGTCACGTCACGTCCGTCCGTCA",
            "AACGTCACGTGGGGCACGTCACGTCACGTCCGTCCGTCA",
            "==========NNNN=========================",
            "10=4N25=", 0, false, false, 0, false);
    logContains("LeftOverlap", "-4..0", 0, 0, 0, 0, 0);
    logContains("LeftSmall", "0..3", 0, 0, 0, 0);
    logContains("LeftLarge", "4..8", 0, 0, 0, 0, 0);
    logContains("RightOverlap", "-4..0", 0, 0, 0, 0, 1);
    logContains("RightSmall", "0..3", 1, 0, 0, 0);
    logContains("RightLarge", "4..8", 1, 0, 0, 0, 0);
  }

  // overlap of 5 (cg_test_errors priors allow max of 4, so it chooses overlap of 4 plus an insert)
  // TODO: correct the cigar generation here (tricky because the overlap region contains a delete!)
  public void testCGOverlap5() {
    checkCG("      acgtc acgtc aacgtcacgtcacgt     acgtcacgtc",
            "ttgca acgtc       aacgtcacgtcacgtgggggacgtcacgtc",
            "ACGTCCAACGTCACGTCACGT.....ACGTCACGTC",
            "ACGT-CAACGTCACGTCACGTGGGGGACGTCACGTC",
            "====IBBBB====================NNNNN==========",
            "4=1I16=5N10=", 2, false, true);
  }

  public void testCGSmallGap2() {
    checkCG("acact ctgggggaaa    aacccccttt       acgtcacgtc",
            "acact   gggggaaa tt aacccccttt ggggg acgtcacgtc",
            "ACACTGGGGGAAA..AACCCCCTTT.....ACGTCACGTC",
            "ACACTGGGGGAAATTAACCCCCTTTGGGGGACGTCACGTC",
            "=====BB==========NN==========NNNNN==========",
            "13=2N10=5N10=", 0, false, true);
  }

  public void testCGSmallGap3() {
    checkCG("acact ctgggggaaa     aacccccttt       acgtcacgtc",
      "acact   gggggaaa ttt aacccccttt ggggg acgtcacgtc",
      "ACACTGGGGGAAA...AACCCCCTTT.....ACGTCACGTC",
      "ACACTGGGGGAAATTTAACCCCCTTTGGGGGACGTCACGTC",
      "=====BB==========NNN==========NNNNN==========",
      "13=3N10=5N10=", 0, false, true);
  }

  public void testCGSmallGap4() {
    // small gap of 4 is too big, so we get a gap of 3 plus 1 delete.
    checkCG("acact ctgggggaaa      aacccccttt       acgtcacgtc",
            "acact   gggggaaa tttt aacccccttt ggggg acgtcacgtc",
            "ACACTGGGGGAAA-...AACCCCCTTT.....ACGTCACGTC".replaceAll("-", ""),
            "ACACTGGGGGAAATTTTAACCCCCTTTGGGGGACGTCACGTC",
            "=====BB==========DNNN==========NNNNN==========",
            "13=1D3N10=5N10=", 2, false, true);
    logContains("LeftOverlap", "-4..0", 0, 0, 1, 0, 0);
    logContains("LeftSmall", "0..3", 0, 0, 0, 1);
    logContains("LeftLarge", "4..8", 0, 1, 0, 0, 0);
    logContains("RightOverlap", "-4..0", 0, 0, 0, 0, 0);
    logContains("RightSmall", "0..3", 0, 0, 0, 0);
    logContains("RightLarge", "4..8", 0, 0, 0, 0, 0);
  }

  public void testCGLargeGap4() {
    checkCG("acact ctgggggaaa aacccccttt      acgtcacgtc",
            "acact   gggggaaa aacccccttt gggg acgtcacgtc",
            "ACACTGGGGGAAAAACCCCCTTT....ACGTCACGTC",
            "ACACTGGGGGAAAAACCCCCTTTGGGGACGTCACGTC",
            "=====BB====================NNNN==========",
            "23=4N10=", 0, false, true);
    assertEquals(1, mSamRecord.getStartPosition());
  }

  public void testCGLargeGap7() {
    checkCG("acact ctgggggaaa aacccccttt         acgtcacgtc",
      "acact   gggggaaa aacccccttt ggggggg acgtcacgtc",
      "ACACTGGGGGAAAAACCCCCTTT.......ACGTCACGTC",
      "ACACTGGGGGAAAAACCCCCTTTGGGGGGGACGTCACGTC",
      "=====BB====================NNNNNNN==========",
      "23=7N10=", 0, false, true);
    assertEquals(1, mSamRecord.getStartPosition());
  }

  public void testCGLargeGap8() {
    checkCG("acact cactgggaaa    aacccccttt         acgtcacgtc",
            "acact     gggaaa ttt aacccccttt gggggggg acgtcacgtc",
            "ACACTGGGAAA...AACCCCCTTT........ACGTCACGTC",
            "ACACTGGGAAATTTAACCCCCTTTGGGGGGGGACGTCACGTC",
            "=====BBBB==========NNN==========NNNNNNNN==========",
            "11=3N10=8N10=", 0, false, true);
    logContains("LeftOverlap", "-4..0", 1, 0, 0, 0, 0);
    logContains("LeftSmall", "0..3", 0, 0, 0, 1);
    logContains("LeftLarge", "4..8", 0, 0, 0, 0, 1);
    logContains("RightOverlap", "-4..0", 0, 0, 0, 0, 0);
    logContains("RightSmall", "0..3", 0, 0, 0, 0);
    logContains("RightLarge", "4..8", 0, 0, 0, 0, 0);
  }

  public void testCGLargeGap9() {
    // large gap of 9 is too big, so we get a gap of 8 plus 1 delete.
    checkCG("acact ctgggggaaa aacccccttt           acgtcacgtc",
            "acact   gggggaaa aacccccttt ggggggggg acgtcacgtc",
            "ACACTGGGGGAAAAACCCCCTTT-........ACGTCACGTC".replaceAll("-", ""),
            "ACACTGGGGGAAAAACCCCCTTTGGGGGGGGGACGTCACGTC",
            "=====BB====================DNNNNNNNN==========",
            "23=1D8N10=", 2, false, true);
    assertEquals(1, mSamRecord.getStartPosition());
  }

  public void testCGComplex1() {
    // several MNPs and indels
    checkCG("acact ctgggACaaa  aac cctttgg        acgtcacgtc",
      "acact   gggggaaa aaacccctttgg gggggg acgtc cgtc",
      "ACACTGGGACAAA..AACCCTTTGG......ACGTCACGTC",
      "ACACTGGGGGAAAAAACCCCTTTGGGGGGGGACGTC-CGTC",
      "=====BB=====XX===NN=X========NNNNNN=====I====",
      "8=2X3=2N1=1X8=6N5=1I4=", 5, false, true);
  }

  public void testCGSoftClippingStart() {
    checkCG("tcaaa aaccggtttt gggggacgtc         cacgtcaagc",
            "  aaa   ccggtttt gggggacgtc actgctt cacgtcaagc",
            "TCAAACCGGTTTTGGGGGACGTC.......CACGTCAAGC",
            "NNAAACCGGTTTTGGGGGACGTCACTGCTTCACGTCAAGC",
            "XX===BB====================NNNNNNN==========",
            "2S21=7N10=", 2, false, true);
    assertEquals(-2, mAlignment.getStart());
  }

  public void testCGSoftClippingEnd() {
    checkCG("tcaaa aaccggtttt gggggacgtc         cacgtcaagc",
            "tcaaa   ccggtttt gggggacgtc actgctt cacgtca",
            "TCAAACCGGTTTTGGGGGACGTC.......CACGTCAAGC",
            "TCAAACCGGTTTTGGGGGACGTCACTGCTTCACGTCANNN",
            "=====BB====================NNNNNNN=======XXX",
            "23=7N7=3S", 3, false, true);
    assertEquals(1, mSamRecord.getStartPosition());
  }

  public void testCGSoftClippingBoth() {
    checkCG("tcaaa aaccggtttt gggggacgtc         cacgtcaagc",
            "  aaa   ccggtttt gggggacgtc actgctt cacgtca",
            "TCAAACCGGTTTTGGGGGACGTC.......CACGTCAAGC",
            "NNAAACCGGTTTTGGGGGACGTCACTGCTTCACGTCANNN",
            "XX===BB====================NNNNNNN=======XXX",
            "2S21=7N7=3S", 5, false, true);
    assertEquals(-2, mAlignment.getStart());
  }

  public void testCgInsertBeforeSpace() {
    checkCG("tatgagaaattctgggttgaaaatt     tttaagaatg",
          "gatatga  aattctgggttgaaaattcttttaagaatgttgaatattggcccccacggtcttctggcttgtagggtttctgcagag",
            "TATGAAATTCTGGGTTGAAAATT....TTTAAGAATG",
            "TATGAAATTCTGGGTTGAAAA--TTCTTTTAAGAATG",
            "=====BB==================IINNNN==========",
            "21=2I4N10=", 3, false, true);
    assertEquals(3, mSamRecord.getStartPosition());
  }

  public void testCgOverlapRealWorld6() {
    checkCG("attcttactanccccaactgagccc     agtaggagta",
            "atactcctacttttgctgggctcagttgggggaagtagaat",
            "ATTCTACTANCCCCAACTGAGCCC......AGTAGGAGTA",
            "ATTCTACTTCCCCCAACTGAGCCCAGCAAAAGTAGGAGTA",
            "=====B====XR==============NNNNNN==========",
            "10=6N14=2X8=", 1, true, true, 0, true, 0, false);
    assertEquals(2, mSamRecord.getStartPosition());
  }

  public void testCgOverlapRealWorld6NsMismatches() {
    checkCG("attcttactanccccaactgagccc     agtaggagta",
            "atactcctacttttgctgggctcagttgggggaagtagaat",
            "ATTCTACTANCCCCAACTGAGCCC......AGTAGGAGTA",
            "ATTCTACTTCCCCCAACTGAGCCCAGCAAAAGTAGGAGTA",
            "=====B====XR==============NNNNNN==========",
            "10=6N14=2X8=", 2, true, true, 0, true, 1, false);
    assertEquals(2, mSamRecord.getStartPosition());
  }

  public void testNgsCgExample() {
    checkCG("TGCATGCATGCATGCACAAAGCTAC      ATTGCTATGC",
            "TGCATGCATGCATGCACATTGCTACAAAAAAATTGCTATGCCCATTACG",
            "TGCATGCATGCATGCACAAAGCTAC......ATTGCTATGC",
            "TGCATGCATGCATGCACATTGCTACAAAAAAATTGCTATGC",
            "==================XX=====NNNNNN==========",
            "18=2X5=6N10=", 2, true, false, 2, true);
  }


  public void testNgsCg2Example() {
    checkCG("TGCATGCATG"
            + "    CATGCATGCACAAAGCTAC",
            "TGCATGCATGCATGCACATTGCTAC",
            "TGCATGCATGCATGCACAAAGCTAC",
            "TGCATGCATGCATGCACATTGCTAC",
            "==========BBBB"
            +     "============XX=====",
            "18=2X5=", 2, true, false, 2, true, 1, true);
  }

  public void testNgsCg2Example2() {
    checkCG(
      "GATCTCTTCACCTCATGATCTACCTGCCT",
//            "AGGCAGGTAGATCATGAGG"
//                           + "TGAAGAGATC",
//            "AGGCAGGTAGATCATGAGGTGAAGAGATC",
      "GATCTCTTCACCTCATGATCTACCTGCCT",
            "GATCTCTTCATCATGATCTACCTGCCT",
            "",
            "==========DDBB"
            + "===================",
            "10=2D17=", 3, true, false, 0,  true, 1, true);
  }

  public void testValidation() {
    checkCG("GTACTGTGTGCGTCCTTGACATCTA     ACCGGCGCCT",
            "GTAGT  GTGCGTCCTTGACATCTAAGGATATCGGCGCCTGA",
            "GTACTGTGCGTCCTTGACATCTA.....ACCGGCGCCT",
            "GTAGTGTGCGTCCTTGACATCTAAGGATATCGGCGCCT",
            "===X=BB====================NNNNN=X========",
            "3=1X19=5N1=1X8=", 2, true, false, 0, true);
  }


  /* This test is dependent on the cg errors priors :/
   public void testBLah() {
    checkCG("TCCATTAAAC.....TTCTTTATAAATTATCCTTCTTTTC",
            "ATGAGTCCATTAAACCTTTTATTCTTTATAAATTATCCAGTCACAGATATTTCTTC",
                 "tccattaaac......ttctttataaattatcct-ttttc",
                 "tccattaaaccttttattctttataaattatccag-tcac",
                 "||||||||||      |||||||||||||||||   |  |",
            "10=6N17=1X1X1=2X1=", 5, false, false, 5, false);

  }*/

  public void testOverlapWeirdScoring() {
    checkCG("tcaaa aagcggtttt gggggacgtc         cacgtcaagc",
            "tcatt ccggtttt gggggacgtc actgctt cacgtcaagcatatat",
            "TCAAAGCGGTTTTGGGGGACGTC.......CACGTCAAGC",
            "TCATTCCGGTTTTGGGGGACGTCACTGCTTCACGTCAAGC",
            "===XXBBXXX=================NNNNNNN==========",
            "3=3X17=7N10=", 5, false, true);
    assertEquals(1, mSamRecord.getStartPosition());    //remember, -4 is 1-based
  }

  public void testLargeSoftClipLeft() {
    checkCG("tcaaa aaccggtttt gggggacgtc         cacgtcaagc",
            "        ccggtttt gggggacgtc actgctt cacgtcaagcatatat",
            "TCAAACCGGTTTTGGGGGACGTC.......CACGTCAAGC",
            "NNNNNCCGGTTTTGGGGGACGTCACTGCTTCACGTCAAGC",
            "XXXXXBBXX==================NNNNNNN==========",
            "5S18=7N10=", 7, false, true);
    assertEquals(-5, mAlignment.getStart());
  }
  public void testLargeSoftClipRight() {
    checkCG("aaaaa aaccggtttt gggggacgtc         cacgtcaagc",
            "aaaaa aaccggtttt gggggacgtc actgctt cacgt",
            "AAAAACCGGTTTTGGGGGACGTC.......CACGTCAAGC",
            "AAAAACCGGTTTTGGGGGACGTCACTGCTTCACGTNNNNN",
            "=====BB====================NNNNNNN=====XXXXX",
        "23=7N5=5S", 5, false, true);
    assertEquals(3, mSamRecord.getStartPosition());    //remember, -4 is 1-based
  }
  public void testLargeSoftClipRight2() {
    checkCG("      cgaactgcac         ctgcaggggg ttttggccaa aaact",
            "tatatacgaactgcac ttcgtca ctgcaggggg ttttggcc",
            "CGAACTGCAC.......CTGCAGGGGGTTTTGGCCAAACT",
            "CGAACTGCACTTCGTCACTGCAGGGGGTTTTGGCCNNNNN",
            "==========NNNNNNN==================XXBBXXXXX",
            "10=7N18=5S", 7, false, false, 7, false);
    assertEquals(7, mSamRecord.getStartPosition());
  }

  public void testEvenLargerSoftLeftClip() {
    checkCG("tcaaa aaccggtttt gggggacgtc         cacgtcaagc",
            "          ggtttt gggggacgtc actgctt cacgtcaagcatatat",
            "TCAAACCGGTTTTGGGGGACGTC.......CACGTCAAGC",
            "NNNNNNNGGTTTTGGGGGACGTCACTGCTTCACGTCAAGC",
            "XXXXXBBXXXX================NNNNNNN==========",
            "7S16=7N10=", 9, true, false, -8, true);
    assertEquals(-7, mAlignment.getStart());
  }

  public void testCgOverlapRealWorld7() {
      checkCG("tcaanaacanggaacagagcagggg     acttgcagaa",
              "ttctgcaagtgttgggcccctgctctgttccctgtttgatatactggggg",
              "TCAANCANGGAACAGAGCAGGGG......ACTTGCAGAA",
              "TCAAACAGGGAACAGAGCAGGGGCCCAACACTTGCAGAA",
              "====RBB====R===============NNNNNN==========",
              "10=6N15=1X2=1X4=", 0, true, true);
  }

  static final String TEMPLATE = ""
    + "GGTGGTCGAGTTGTGTGTAGGCGATGGCACGACCTCAGGAGTTACTGTACCATGGTTGCTTCTAGCGAGAAATTCTTG"
    ;
  static final String[] READS = {
    //One-Based-Start,  Read-String,                  Cigar
    "6",  "TCGAGAGTTGTGTGTAGGCGATGGC     TCAGGAGTTA".replaceAll(" ", ""), "23=6N10=",
    "16", "TGTAGAGGCGATGGCACCTCAGGAG     ACCATGGTTG".replaceAll(" ", ""), "13=3N10=7N10=",
    "21", "GCGATGATGGCACGACCTCAAGAGT     ACCATGGTTG".replaceAll(" ", ""), "17=1X4=6N10=",
    "22", "CGATGTGGCACGACCTCAGGAGTTA     CATGGTTGCT".replaceAll(" ", ""), "23=6N10=",
    "40", "AGTTATTCTGTACCATGGTTGCTTC     AATTCTTGCA".replaceAll(" ", ""), "23=8N8=2S"
  };
  public void testMultiReads() throws InvalidParamsException, IOException {
    final ByteArrayOutputStream bos = new ByteArrayOutputStream();
    Diagnostic.setLogStream(new PrintStream(bos));
    final RealignParamsImplementation params = new RealignParamsImplementation(MachineErrorParams.builder().errors("cg_test_errors").create());
    final CgGotohEditDistance ed = new CgGotohEditDistance(7, params, 0);
    final byte[] tmpl = DnaUtils.encodeString(TEMPLATE);
    for (int i = 0; i < READS.length; i += 3) {
      final int start = Integer.parseInt(READS[i]);
      final byte[] read = DnaUtils.encodeString(READS[i + 1]);
      final String cigar = READS[i + 2];
      final int approxStart = start ^ 5; // move the start around a bit
      final int[] v = ed.calculateEditDistance(read, read.length, tmpl, approxStart, Integer.MAX_VALUE, 7, true);
      //      System.out.println(ed.toString());
      //      System.out.println("start: " + ActionsHelper.zeroBasedTemplateStart(v)
      //          + " score: " + ActionsHelper.alignmentScore(v));
      //      System.out.println("actions: " + ActionsHelper.toString(v));
      final AlignmentResult alignment = new AlignmentResult(read, v, tmpl);
      alignment.setIdentifyingInfo(true, false);
      assertEquals(cigar, alignment.getCigarString(false, false));
      assertEquals(start - 1, ActionsHelper.zeroBasedTemplateStart(v));
      assertEquals(start - 1, alignment.getStart());
    }
    ed.logStats();
    mLog = bos.toString();
    logContains("LeftOverlap", "-4..0", 0, 1, 4, 0, 0);
    logContains("LeftSmall", "0..3", 4, 0, 0, 1);
    logContains("LeftLarge", "4..8", 0, 0, 3, 1, 1);
    logContains("RightOverlap", "-4..0", 0, 0, 0, 0, 0);
    logContains("RightSmall", "0..3", 0, 0, 0, 0);
    logContains("RightLarge", "4..8", 0, 0, 0, 0, 0);
    Diagnostic.setLogStream();
  }

  public static void main(String[] args) {
    new CgGotohEditDistanceTest().testCGLeftRightSymmetry();
  }

  public void testCGLeftRightSymmetry() {
    getEditDistanceInstance(0); // new CgGotohEditDistance(7, new ScoreMatrixCGTest.MockRealignParamsCG());
    final CgGotohEditDistance ed = mCgGotoh;
    final String[] reads = {
        "ATAAA AAGGCGACAT GCCAATGTGT        TTCAACTTTC",
        "ACGTT CGTTACGTAC GTACGTACGT        ACGACGACGA",
        "GGGAA ACATATATAT GGCGACATGC        ATTTTTACAC",
        "GGGTG TGCAATTGAC ACGACGACGA        TTCAACTTTC",
        "ATCCA ACGGCCACAT GCCATTTTGT        CGTTACGTAC",
        "AAAAA AAAAAAAAAA AAAAAAAAAA        AAAAAAAAAA",
        "AAAAA AACACACACA CACACACACA        ACACACACAA",
        "AAAAA CCCCCCCCCC GGGGGGGGGG        TTTTTTTTTT",
        "GACGT ACGTATACGT GGTGACTTGC        ATACCCACGC",
        "GACGT GCATACTATC GTATAAAAGC        TCGCCCACTT",
    };
    final PortableRandom rand = new PortableRandom(42);
    final String[] mutations = {"A", "C", "G", "T", "", "AC", "GT"};
    String mut = "";
    for (int repeat = 0; repeat < reads.length; ++repeat) {
      final String read = reads[repeat % reads.length].replaceAll(" ", "").toLowerCase(Locale.ROOT);
      String tmpl = "aaaagggg" + read.substring(0, 25) + "tatata" + read.substring(25) + "accaacca";
      //System.out.println("@@@@@@@@@@@@@@@ starting with read:" + read + "@@@@@@@@@@@@@@@@@@@");
      int readStart = 8;
      int readEnd = readStart + 39;  // ie. 35 - 2 + 0 + 6 (for the most common gap sizes)
      for (int age = 0; age < 10 && Math.abs(readEnd - readStart - 39) < 5 ; ++age) {
        final byte[] s1 = DnaUtils.encodeString(read);
        final byte[] s2 = DnaUtils.encodeString(tmpl.toLowerCase(Locale.ROOT));
        final int[] tmp = ed.calculateEditDistance(s1, s1.length, s2, readStart, Integer.MAX_VALUE, 7, true);
        final int[] v = new int[tmp.length];
        System.arraycopy(tmp, 0, v, 0, tmp.length);
        //      System.out.println(ed.toString());
        //      System.out.println("start: " + ActionsHelper.zeroBasedTemplateStart(v)
        //          + " score: " + ActionsHelper.alignmentScore(v));
        //      System.out.println("actions: " + ActionsHelper.toString(v));
        final int leftScore = ActionsHelper.alignmentScore(v);

        // now check that we get the same alignment score if we flip read+template onto the opposite arm.
        final byte[] s1R = DnaUtils.encodeString(StringUtils.reverse(read).replaceAll(" ", "").toLowerCase(Locale.ROOT));
        final byte[] s2R = DnaUtils.encodeString(StringUtils.reverse(tmpl).replaceAll(" ", "").toLowerCase(Locale.ROOT));
        assertEquals(s1.length, s1R.length);
        assertEquals(s2.length, s2R.length);
        final int startR = s2.length - readEnd ;
        final int[] vR = ed.calculateEditDistance(s1R, s1R.length, s2R, startR, Integer.MAX_VALUE, 7, false);
        final int rightScore = ActionsHelper.alignmentScore(vR);
        //System.out.println(mCgGotoh.toString());
        //assertTrue(Math.abs(leftScore - rightScore) <= 1);  // but they should be equal!
        final int diff = leftScore - rightScore;
        final String actsLeft = ActionsHelper.toString(v);
        final String actsRight = ActionsHelper.toString(vR);
        final StringBuilder sb = new StringBuilder(actsRight);
        sb.reverse();
        if (diff != 0) {
          System.out.println("read=" + StringUtils.spaces(readStart + 1) + read);
          System.out.println("tmpl=" + tmpl.substring(0, readStart) + " " + tmpl.substring(readStart, readEnd) + " " + tmpl.substring(readEnd));
          System.out.println("" + actsLeft + " " + leftScore + " " + rightScore + " " + sb.toString() + "  last mutation: " + mut + " age=" + age);
        }
        assertEquals(leftScore, rightScore);
        assertEquals(actsLeft, sb.toString());
        // now mutate the template a bit
        final int r = rand.nextInt(mutations.length * (tmpl.length() - 2));
        mut = mutations[r % mutations.length];
        final int pos = r / mutations.length;
        if (pos <= readStart) {
          readStart = readStart - 1 + mut.length();
        }
        if (pos <= readEnd) {
          readEnd = readEnd - 1 + mut.length();
        }
        tmpl = tmpl.substring(0, pos) + mut + tmpl.substring(pos + 1);
      }
    }
    //System.out.println("Score Differences: " + java.util.Arrays.toString(hist));
  }

  /* TODO: add in some of these tests too.  (They usually need some modification).
  public void testOverlap2Ugly() {

      checkCG("tatccgaaattctgggttgaaaatt     tttaagaatg",
              "gatatga  aattctgggttgaaaatt  cttttaagaatgttgaatattggcccccacggtcttctggcttgtagggtttctgcagag",
              "tatccgaaattctgggttgaaaatt.....tttaagaatg",
              "tat-.gaaattctgggttgaaaa--ttctt-ttaagaatg",
              "|||  ||||||||||||||||||        |||||||||",
              "3=1I18=2I5N1I9=", 7, false, true);
  }
  public void testOverlap2() {
      check("tatgagaaattctgggttgaaaatt     tttaagaatg",
            "gatatgaaattctgggttgaaaattcttttaagaatgttgaatattggcccccacggtcttctggcttgtagggtttctgcagag",
            "tatgagaaattctgggttgaaaatt.....tttaagaatg",
            "tat..gaaattctgggttgaaaa--ttctt-ttaagaatg",
            "|||  ||||||||||||||||||        |||||||||",
            "21=2I5N1I9=", 5, false);
  }
  public void testOverlap2RC() {
      check("tatgagaaattctgggttgaaaatt     tttaagaatg",
            "ccaatattcaacattcttaaaagaattttcaacccagaatttcatatc",
            "tatgagaaattctgggttgaaaatt.....tttaagaatg",
            "tat..gaaattctgggttgaaaa--ttctt-ttaagaatg",
            "|||  ||||||||||||||||||        |||||||||",
            "9=1I5N2I21=", 5, true);
  }

  public void testOverlap2Rev() {
      check("gtaagaattt     ttaaaagttgggtcttaaagagtat",
            "gtaagaattttcttaaaagttgggtcttaaagtatag",
            "gtaagaattt.....ttaaaagttgggtcttaaagagtat",
            "gtaagaa---ttttcttaaaagttgggtcttaaag..tat",
            "|||||||        ||||||||||||||||||||  |||",
            "7=3I5N23=", 4, false);
  }

  public void testCgOverlap1() {
      check("tatcaaaaattctgggttgaaaatt     tttaagaatg",
            "gatatcaaaattctgggttgaaaattcttttaagaatgttg",
            "tatcaaaaattctgggttgaaaatt.....tttaagaatg",
            "tatc.aaaattctgggttgaaaa--ttctt-ttaagaatg",
            "|||| ||||||||||||||||||        |||||||||",
            "22=2I5N1I9=", 5, false);
  }

  public void testCgOverlap1Rev() {
      check("gtaagaattt     ttaaaagttgggtcttaaaaactat",
            "gttgtaagaattttcttaaaagttgggtcttaaaactatag",
            "gtaagaattt.....ttaaaagttgggtcttaaaaactat",
            "gtaagaa---ttttcttaaaagttgggtcttaaaa.ctat",
            "|||||||        |||||||||||||||||||| ||||",
            "7=3I5N24=", 4, false);
  }
  public void testCgOverlap1RevRC() {
      check("gtaagaattt     ttaaaagttgggtcttaaaaactat",
            "ctatagttttaagacccaacttttaagaaaattcttacaac",
            "gtaagaattt.....ttaaaagttgggtcttaaaaactat",
            "gtaagaa---ttttcttaaaagttgggtcttaaaa.ctat",
            "|||||||        |||||||||||||||||||| ||||",
            "24=5N3I7=", 4, true);
  }

  //  not supporting overlap of 3 properly yet.
  public void testCgOverlap3() {
    final AlignmentResult ar =
      check("atgaagaaattctgggttgaaaatt     tttaagaatg",
            "gatatgaaattctgggttgaaaattcttttaagaatgttgaatat",
            "atgaagaaattctgggttgaaaatt.....tttaagaatg",
            "at---gaaattctgggttgaaaa--ttctt-ttaagaatg",
            "||   ||||||||||||||||||        |||||||||", 5, false, 3);
    final String cigar = CigarFormatter.alignmentToCigar(ar.mRead, ar.mTemplate, TOTAL_REFLENGTH_LONG, ar.mMatch, false, ar.mStart);
    assertEquals("2M3I18M2I5N1I9M", cigar);
  }

  // Test that '.' in a read is treated the same as 'n'.
  public void testAlignWithIndelsCG5() {
    final String read = "tccttag     gatgcaagacaagagggcctctc";
    final byte[] rr = DnaUtils.encodeArray(read.getBytes());
    EditDistance ed = getEditDistanceInstance(1, 1);
    final String t = "aaaaaatccttagactcagaggatgcaagacaagagggcctctcgaatgt";
                          //tccttag   .....gatgcaagacaagagggcctctc
    final byte[] template = DnaUtils.encodeArray(t.getBytes());
    int[] actions = ed.calculateEditDistance(rr, rr.length, template, 7, false, Integer.MAX_VALUE, 9);
    final AlignmentResult ar = new AlignmentResult(rr, actions, template.length, template, false);

    ar.setIdentifyingInfo(true, false);
    ar.setRemainingOutput(1, "toxic");
    final SAMRecord sr = ar.toSamRecord(header(), false, null);
    final String srformatted = sr.format();
    assertTrue(srformatted.contains("TCCTTAGGATGCAAGACAAGAGGGCCTCTC"));
    assertTrue(srformatted.contains("7=7N1D23="));
    assertTrue(srformatted.contains("AS:i:2\tNM:i:1"));

    final String rcTemplate = DnaUtils.reverseComplement(t);
    ed = getEditDistanceInstance(1, 1);
    final byte[] rctemplateBytes = DnaUtils.encodeArray(rcTemplate.getBytes());
    actions = ed.calculateEditDistance(rr, rr.length, rctemplateBytes, 7, true, Integer.MAX_VALUE, 9);
    final AlignmentResult ar2 = new AlignmentResult(rr, actions, rctemplateBytes.length, rctemplateBytes, false);

    ar2.setIdentifyingInfo(true, true);
    ar2.setRemainingOutput(1, "toxic");
    final SAMRecord sr2 = ar2.toSamRecord(header(), false, null);

    final String sr2formatted = sr2.format();
    assertTrue(sr2formatted.contains("GAGAGGCCCTCTTGTCTTGCATCCTAAGGA"));
    assertTrue(sr2formatted.contains("23=1D7N7="));
    assertTrue(sr2formatted.contains("AS:i:2\tNM:i:1"));
  }

  public void testCgOverlapRealWorld() {
      check("ataaaaaggcgacatgccaatgtgt     ttcaactttc",
            "ataaaaaaggcgacatgccaatgtgtcgcctttttcaactttccgattaagaacctgctcagcgggt",
            "ataaaaaggcgacatgccaatgtgt.......ttcaactttc",
            "aaa..aaggcgacatgccaatgtgtcgcctttttcaactttc",
            "| |  ||||||||||||||||||||       ||||||||||",
            "1=1X21=7N10=", 1, false);
  }
  public void testCgOverlapRealWorld2() {
      check("taacgccgat     ccaacggttatctcgatttttttta",
            "taacgccgattgaggccaacggttatctcgatttttttatc",
            "taacgccgat.....ccaacggttatctcgatttttttta",
            "taacgccgattgaggccaacggttatctcgatttt.ttta",
            "||||||||||     |||||||||||||||||||| ||||",
            "10=5N24=", 0, false);
  }

  public void testCgOverlapRealWorld3() {
      check("cgacggcggt     tgcttcactatgggaaagaggtcca",
            "cgacggcggtgggattgcttcactatgggaaagagtcacat",
            "cgacggcggt.....tgcttcactatgggaaagaggtc-ca",
            "cgacggcggtgggattgcttcactatgggaaagag.tcaca",
            "||||||||||     |||||||||||||||||||| || ||",
            "10=5N22=1D2=", 2, false);
  }

  public void testCgOverlapRealWorld4() {
      check("tttgtgtaggtcggataaggcgttc     atccgacacg",
            "ttttgtaggtcggataaggcgttcacgccgcatccgacacgg",
            "tttgtgtaggtcggataaggcgttc.......atccgacacg",
            "ttt..gtaggtcggataaggcgttcacgccgcatccgacacg",
            "|||  ||||||||||||||||||||       ||||||||||",
            "23=7N10=", 0, false, 1);
  }

  public void testCgOverlapRealWorld5() {
    final AlignmentResult ar =
      check("cgacggcggt     tgcttcactatgggaaagaggtcca",
            "cgacggcggtgggattgcttcactatgggaaagagtcacat",
            "cgacggcggt.....tgcttcactatgggaaagaggtc-ca",
            "cgacggcggtgggattgcttcactatgggaaagag.tcaca",
            "||||||||||     |||||||||||||||||||| || ||",
            "10=5N22=1D2=", 2, false, 1);
    assertEquals(1, ar.mismatches());
  }



  public void testCgRegression() {
    final AlignmentResult ar =
      check("cnggttaaaatatgaagtgaccacc     atgcttgaga",
            "ctggtaaaatatgaagtgaccaccaaagggagcttgagagaggagaaaatgact",
            "cnggttaaaatatgaagtgaccacc.....atgcttgaga",
            "ctgg.taaaatatgaagtgaccaccaaagggagcttgaga",
            "|||| ||||||||||||||||||||       ||||||||",
            "24=5N2X8=", 2, false, 1);
    assertEquals(2, ar.mismatches());
  }

  public void testFromAlignmentResult() {
    final AlignmentResult ar =
    check("tagacaaatg     ggattacaagccacaggagggggaa",
          "tagacaaatgtgactggattacaagccacaggaga",
          "tagacaaatg.....ggattacaagccacaggagggggaa",
          "tagacaaatgtgactggattacaagccacaggaga..nnn",
          "||||||||||     |||||||||||||||||||      ",
            "10=5N19=1X3S", 4, false, 1);
    assertEquals(4, ar.mismatches());
  }

  public void testFromAlignmentResult2() {
    final AlignmentResult ar =
    check("ctgctgaccc     gacctactggcggactcgggggtta",
            "ctgctgacccagagccggacctactggcggactcgggtta",
            "ctgctgaccc.......gacctactggcggactcgggggtta",
            "ctgctgacccagagccggacctactggcggactcggg..tta",
            "||||||||||       ||||||||||||||||||||  |||",
            "10=7N23=", 0, false, 1);
    assertEquals(0, ar.mismatches());
  }

  public void testBug609SoftClippingEnd() {
    final AlignmentResult ar =
      check("tagacaaatg     ggattacaagccacaggagtgaata",
            "tagacaaatgtgactggattacaagccacaggaga",
            "tagacaaatg.....ggattacaagccacaggagtgaata",
            "tagacaaatgtgactggattacaagccacaggaga.nnnn",
            "||||||||||     |||||||||||||||||||      ",
            "10=5N19=1X4S", 5, false);
    assertEquals(5, ar.mismatches());
  }

  public void testBug609SoftClippingStart() {
    final AlignmentResult ar =
      check("ttctgtggggtgactggattacaag     ggagtgaata",
            "tggggtgactggattacaagaaaaaggagtgaat",
            "ttctgtggggtgactggattacaag.....ggagtgaata",
            "nnn..tggggtgactggattacaagaaaaaggagtgaatn",
            "     ||||||||||||||||||||     ||||||||| ",
            "3S20=5N9=1X", 4, false, -3);
    assertEquals(4, ar.mismatches());
  }

  private AlignmentResult getAlignment(final String read, final String template, final EditDistance f, final boolean rc) {
    return getAlignment(read, template, f, 0, rc);
  }

  private AlignmentResult getAlignment(final String read, final String template, final EditDistance f, final int startPos, final boolean rc) {
    final byte[] s1 = read.getBytes();
    final byte[] s2 = template.getBytes();
    DnaUtils.encodeArray(s1);
    DnaUtils.encodeArray(s2);
    final int[] actions = f.calculateEditDistance(s1, s1.length, s2, startPos, rc, Integer.MAX_VALUE, 9);
    final AlignmentResult alignment = new AlignmentResult(s1, actions, s2.length, s2, false);
    if (rc) {
      alignment.setIdentifyingInfo(alignment.isFirst(), true);
    }
    return alignment;
  }
  public void testMultipleCgAlignments() {
    // Check that any object churn reduction isn't overriding something.
    final EditDistance f = getEditDistanceInstance(1, 1);
    final AlignmentResult alignment3 = getAlignment("attcttactanccccaactgagccc     agtaggagta",
            "atactcctacttttgctgggctcagttgggggaagtagaat", f, true);
    final String read1 = "cgacggcggt.....tgcttcactatgggaaagaggtcca";
    final AlignmentResult alignment1 = getAlignment(read1, "cgacggcggtgggattgcttcactatgggaaagagtcacat", f, false);
    final String template = "ctatagttttaagacccaacttttaagaaaattcttacaac";
    final String read2 = "gtaagaattt.....ttaaaagttgggtcttaaaaactat";
    final AlignmentResult alignment2 = getAlignment(read2, template, f, true);
    alignment2.setIdentifyingInfo(true, true);
    final String cigar1 = alignment1.getCigarString(alignment1.getStart(), false);
    final String cigar2 = alignment2.getCigarString(alignment2.getStart(), true);
    assertEquals("10=5N22=1D2=", cigar1);
    assertEquals("24=5N3I7=", cigar2);
    assertEquals("gtaagaattt.....ttaaaagttgggtcttaaaaactat" + "\t"
            + "gtaagaa---ttttcttaaaagttgggtcttaaaa.ctat" + "\t"
            + "|||||||        |||||||||||||||||||| ||||" , alignment2.tabularString());
    assertEquals("cgacggcggt.....tgcttcactatgggaaagaggtc-ca" + "\t"
            + "cgacggcggtgggattgcttcactatgggaaagag.tcaca" + "\t"
            + "||||||||||     |||||||||||||||||||| || ||", alignment1.tabularString());
    assertEquals("attcttactanccccaactgagccc......agtaggagta" + "\t"
            + "attc.tacttcccccaactgagcccagcaaaagtaggagta" + "\t"
            + "|||| |||| |||||||||||||||      ||||||||||", alignment3.tabularString());
  }
   */

  /* public void testRcBullcrap() {

    final EditDistance ed = new CgEditDistance(1);
    final AlignmentResult ar2 = getAlignment("attaaaaaaattttttttttttttt.....gacaaagtct", "cgacagagcgagattccgtctcaaacaaaaaaaaaaaaaaaaaaaattagctgggtgtgattataggtgc", ed, 16, true);
      //                        gcacctataatcacacccagctaattttttttttttttttttttgtttgagacggaatctcgctctgtcg < 14 start                           x
      //                   gcacctataatcacacccagctaattttttttttttttttttttgtttgagacggaatctcgctctgtcg < 19 start
    System.err.println(ar2.getScore());
    System.err.println(ar2.getActionsString());

    System.err.println("-------------------------------------------------");
    System.err.println("-------------------------------------------------");
    System.err.println("-------------------------------------------------");
    System.err.println("-------------------------------------------------");

    final EditDistance ed2 = new LoopingEditDistance(new RcEditDistance(new CgEditDistance(1)));
    final AlignmentResult ar3 = getAlignment("attaaaaaaattttttttttttttt.....gacaaagtct", "cgacagagcgagattccgtctcaaacaaaaaaaaaaaaaaaaaaaattagctgggtgtgattataggtgc", ed2, 16, true);
    System.err.println(ar3.getScore());
    System.err.println(ar3.getActionsString());
    assertEquals(8, ar3.getScore());

    System.err.println("-------------------------------------------------");
    System.err.println("-------------------------------------------------");
    System.err.println("-------------------------------------------------");
    System.err.println("-------------------------------------------------");


    final EditDistance ed3 = new RcEditDistance(new LoopingEditDistance(new CgEditDistance(1)));
    final AlignmentResult ar4 = getAlignment("attaaaaaaattttttttttttttt.....gacaaagtct", "cgacagagcgagattccgtctcaaacaaaaaaaaaaaaaaaaaaaattagctgggtgtgattataggtgc", ed3, 16, true);
    System.err.println(ar4.getScore());
    System.err.println(ar4.getActionsString());

    assertEquals(8, ar4.getScore());
  }*/

  public void testRowOffsetLeftArm() {
    getEditDistanceInstance(0);
    final int[] offsets = mCgGotoh.makeRowOffsets(CgGotohEditDistance.LEFT_ARM);
    assertEquals(36, offsets.length);
    final int half = -8; // half of the width of the band (at row 1).
    assertEquals(half, offsets[0]);
    assertEquals(half + 1, offsets[1]);
    assertEquals(half + 5, offsets[5]);
    assertEquals(half + 6 - 2, offsets[6]);
    assertEquals(half + 15 - 2, offsets[15]);
    assertEquals(half + 16 - 2, offsets[16]);
    assertEquals(half + 25 - 2, offsets[25]);
    assertEquals(half + 26 - 2 + 6, offsets[26]);
    assertEquals(half + 35 - 2 + 6 , offsets[35]);
  }

  public void testRowOffsetRightArm() {
    getEditDistanceInstance(0);
    final int[] offsets = mCgGotoh.makeRowOffsets(CgGotohEditDistance.RIGHT_ARM);
    assertEquals(36, offsets.length);
    final int half = -8; // half of the width of the band (at row 1).
    assertEquals(half, offsets[0]);
    assertEquals(half + 1, offsets[1]);
    assertEquals(half + 10, offsets[10]);
    assertEquals(half + 11 + 6, offsets[11]);
    assertEquals(half + 20 + 6, offsets[20]);
    assertEquals(half + 21 + 6, offsets[21]);
    assertEquals(half + 30 + 6, offsets[30]);
    assertEquals(half + 31 + 6 - 2, offsets[31]);
    assertEquals(half + 35 + 6 - 2, offsets[35]);
  }

  public void testFixed() {
    getEditDistanceInstance(0);
    assertNull(mCgGotoh.calculateEditDistanceFixedStart(null, 0, 0, null, 0, 0, 0));
    assertNull(mCgGotoh.calculateEditDistanceFixedEnd(null, 0, 0, null, 0, 0, 0, 0));
    assertNull(mCgGotoh.calculateEditDistanceFixedBoth(null, 0, 0, null, 0, 0, 0, 0));
  }

  public void testSoftClip() {
    final MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());
    try {
      getEditDistanceInstance(0);

      final byte[] t = DnaUtils.encodeString("gatcatcggcgacatgccattgtgttttttttcaacttt");
      final byte[] s1 = DnaUtils.encodeString("gatcatcggcgacatgccattgtgt     ttcaactttc".replaceAll(" ", ""));
      final byte[] t2 = DnaUtils.encodeString("gatcatcggcgacatgccattgtgttttttttcaa");

      //System.err.println(s1.length);
      int[] actions = mCgGotoh.calculateEditDistance(s1, s1.length, t, 0, 10, 5, true);
      if (actions != null) {
        assertEquals("=========================NNNNN=========X", ActionsHelper.toString(actions));
        assertEquals(1, actionsAlignmentScore(actions, false));
        assertEquals(0, actions[ActionsHelper.TEMPLATE_START_INDEX]);
      }

      actions = mCgGotoh.calculateEditDistance(s1, s1.length, t2, 0, 10, 5, true);
      if (actions != null) {
        assertEquals("=========================NNNNN=====XXXXX", ActionsHelper.toString(actions));
        assertEquals(5, actionsAlignmentScore(actions, false));
        assertEquals(0, actions[ActionsHelper.TEMPLATE_START_INDEX]);
      }

      getEditDistanceInstance(1);

      actions = mCgGotoh.calculateEditDistance(s1, s1.length, t, 0, 10, 5, true);
      if (actions != null) {
        assertEquals(1, actionsAlignmentScore(actions, false));
        assertEquals(0, actions[ActionsHelper.TEMPLATE_START_INDEX]);
        assertEquals("=========================NNNNN=========X", ActionsHelper.toString(actions));
      }

      actions = mCgGotoh.calculateEditDistance(s1, s1.length, t2, 0, 10, 5, true);
      if (actions != null) {
        assertEquals(5, actionsAlignmentScore(actions, false));
        assertEquals(0, actions[ActionsHelper.TEMPLATE_START_INDEX]);
        assertEquals("=========================NNNNN=====XXXXX", ActionsHelper.toString(actions));
      }
    } finally {
      Diagnostic.setLogStream();
    }
  }

  public void testalign3Mismatches() {
    getEditDistanceInstance(0);

    final byte[] r = DnaUtils.encodeString("gatcatatatCtatatatatatata     tatataaaac".replaceAll(" ", ""));
    final byte[] t = DnaUtils.encodeString("gatcatataCatatatatatCtataNNNNNtatataaaac");

    int[] actions = mCgGotoh.calculateEditDistance(r, r.length, t, 0, 10, 3, true);
    if (actions != null) {
      assertEquals("=========XX=========X====NNNNN==========", ActionsHelper.toString(actions));
      assertEquals(3, actionsAlignmentScore(actions, false));
      assertEquals(0, actions[ActionsHelper.TEMPLATE_START_INDEX]);
    }
    getEditDistanceInstance(1);
    actions = mCgGotoh.calculateEditDistance(r, r.length, t, 0, 10, 3, true);
    if (actions != null) {
      assertEquals("=========XX=========X====NNNNN==========", ActionsHelper.toString(actions));
      assertEquals(3, actionsAlignmentScore(actions, false));
      assertEquals(0, actions[ActionsHelper.TEMPLATE_START_INDEX]);
    }
  }

  public void testAlignNs() {
    getEditDistanceInstance(0);

    final byte[] r = DnaUtils.encodeString("gatcatatatNtatatatatNtata     tatataaaac".replaceAll(" ", ""));
    final byte[] t = DnaUtils.encodeString("gatcatataNatatatatatNtataNNNNNtatataaaac");

    int[] actions = mCgGotoh.calculateEditDistance(r, r.length, t, 0, 10, 3, true);
    if (actions != null) {
      assertEquals("=========TR=========R====NNNNN==========", ActionsHelper.toString(actions));
      assertEquals(0, actionsAlignmentScore(actions, false));
      assertEquals(0, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
      assertEquals(0, actions[ActionsHelper.TEMPLATE_START_INDEX]);
    }
    getEditDistanceInstance(1);
    actions = mCgGotoh.calculateEditDistance(r, r.length, t, 0, 10, 3, true);
    if (actions != null) {
      assertEquals("=========TR=========R====NNNNN==========", ActionsHelper.toString(actions));
      assertEquals(3, actionsAlignmentScore(actions, true));
      assertEquals(3, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
      assertEquals(0, actions[ActionsHelper.TEMPLATE_START_INDEX]);
    }
  }

  public void testAlignmentScores() {

    getEditDistanceInstance(0, "cg_test_errors", false);
//    Claimed alignment is incorrect:
//      ACCTCCTATGAAAAAACTTCCTACCACTCAC-CCTAG <-- template
//      ||||||||||      ||||||||||  ||| ||| |
//      ACCTCCTATG......CTTCCTACCAAACACTNCTTG <-- read
//      Observed mismatches: 4 Claimed mismatches: null
//      Observed alignment score: 5 Claimed alignment score: 3
//      Alignment mismatch 644988       137     chr1    554579  0       10=6N10=2X3=1I3=1X1=    *       0       0       ACCTCCTATGCTTCCTACCAAACACTNCTTG 6877-4.&04;:::;;7;:::89998!9988 AS:i:3
//    XS:Z:10=6N10=2X3=1I1X2I1=2B3=1X1=       XR:Z:AATGATT    XQ:Z:5578       IH:i:2  NH:i:2

    byte[] r = DnaUtils.encodeString("ACCTCCTATG     CTTCCTACCAAACACTGATCNCTTG".replaceAll(" ", ""));
    byte[] t = DnaUtils.encodeString("ACCTCCTATGAAAAAACTTCCTACCACTCACCCTAGGGGGGGGGGGG");
    int[] actions = mCgGotoh.calculateEditDistance(r, r.length, t, 0, 10, 3, false);

//    System.err.println(ActionsHelper.toString(actions) + " as=" + ActionsHelper.alignmentScore(actions));
    assertEquals(9, ActionsHelper.alignmentScore(actions));




//    GAGGAGGAGGGAATGTTTCCAAACTCAGATTTCTATGAGGCC <-- template
//    || || ||||| |   ||||||||||      ||||||||||
//    GATGNTGAGGGTAA..TTCCAAACTC......CTATGAGGCC <-- read
//    Observed mismatches: 4 Claimed mismatches: null
//    Observed alignment score: 4 Claimed alignment score: 3
//    Alignment mismatch 727966       153     chr22   45958796        37      2=1X2=1X5=1X1=1X2N10=6N10=      *       0       0 GATGNTGAGGGTAATTCCAAACTCCTATGAGGCC       899:!87989999::;;<;:98:9-0*2/74384      AS:i:3
//    XS:Z:2=1X2=2B1=1D1X5=1X1=1X2N10=6N10=      XR:Z:TTTA       XQ:Z:6  IH:i:1  NH:i:1
    r = DnaUtils.encodeString("GATGNGTGAGGGTAATTCCAAACTC     CTATGAGGCC".replaceAll(" ", ""));
    t = DnaUtils.encodeString("GAGGAGGAGGGAATGTTTCCAAACTCAGATTTCTATGAGGCC");

    actions = mCgGotoh.calculateEditDistance(r, r.length, t, 0, 10, 3, true);

//    System.err.println(ActionsHelper.toString(actions) + " as=" + ActionsHelper.alignmentScore(actions));
    assertEquals(6, ActionsHelper.alignmentScore(actions));


//    Claimed alignment is incorrect:
//      GAAGGTTGGGGGAGTGGGAGTTGGTGGCCTACCACTG-GGT <-- template
//      ||||||||||       |||||||||||||||||||| |||
//      GAAGGTTGNG.......GAGTTGGTNGCCTACCACTGTGGT <-- read
//      Observed alignment score: 4 Claimed alignment score: 3
//      Alignment mismatch 62705        131     paolo-bac       698     55      8=1X1=7N8=1X11=1I3=     paolo-bac       363     -335    GAAGGTTGNGGAGTTGGTNGCCTACCACTGTGGT      /21-/704!3.062/510!1,,,1042+/-5021      AS:i:3  NM:i:3  MQ:i:255        GS:Z:GG GC:Z:29S1G4S    GQ:Z:0  XA:i:3  IH:i:1  NH:i:1
    r = DnaUtils.encodeString("GAAGGTTGNG     GAGTTGGTNGCCTACCACTGGTGGT".replaceAll(" ", ""));
    t = DnaUtils.encodeString("GAAGGTTGGGGGAGTGGGAGTTGGTGGCCTACCACTGGGT");


    getEditDistanceInstance(1);
    actions = mCgGotoh.calculateEditDistance(r, r.length, t, 0, 10, 3, false);

//    System.err.println(ActionsHelper.toString(actions) + " as=" + ActionsHelper.alignmentScore(actions));
    assertEquals(4, ActionsHelper.alignmentScore(actions));

  }

  private int actionsAlignmentScore(int[] actions, boolean treatNsAsMismatches) {
    final ActionsHelper.CommandIterator it = ActionsHelper.iterator(actions);
    int prevAction = ActionsHelper.SAME;
    int ascore = 0;
    while (it.hasNext()) {
      final int currAction = it.next();
      if (currAction == ActionsHelper.MISMATCH) {
        ++ascore;
      } else if (currAction == ActionsHelper.DELETION_FROM_REFERENCE) {
        if (prevAction != ActionsHelper.DELETION_FROM_REFERENCE) {
          ascore += 2;
        } else {
          ++ascore;
        }
      } else if (currAction == ActionsHelper.INSERTION_INTO_REFERENCE) {
        if (prevAction != ActionsHelper.INSERTION_INTO_REFERENCE) {
          ascore += 2;
        } else {
          ++ascore;
        }
      } else if (currAction == ActionsHelper.UNKNOWN_TEMPLATE || currAction == ActionsHelper.UNKNOWN_READ) {
        if (treatNsAsMismatches) {
          ++ascore;
        }
      }
      prevAction = currAction;
    }
    return ascore;
  }
}

