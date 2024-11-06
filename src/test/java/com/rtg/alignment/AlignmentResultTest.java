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
package com.rtg.alignment;

import java.io.File;

import com.rtg.mode.DnaUtils;
import com.rtg.ngs.NgsParams;
import com.rtg.ngs.NgsParamsBuilder;
import com.rtg.ngs.tempstage.BinaryTempFileRecord;
import com.rtg.sam.SamBamConstants;
import com.rtg.sam.SamUtils;
import com.rtg.sam.SamValidator;
import com.rtg.util.io.MemoryPrintStream;

import htsjdk.samtools.SAMRecord;
import junit.framework.TestCase;

/**
 */
public class AlignmentResultTest extends TestCase {

  private byte[] getSequence(final String s) {
    return DnaUtils.encodeArray(s.replaceAll(" ", "").getBytes());
  }

  public void testMismatches() {
    AlignmentResult ar = new AlignmentResult(DnaUtils.encodeString("acgtacgta"), ActionsHelper.build("=========", 0, 0), DnaUtils.encodeString("acgtacgta"));
    assertEquals(0, ar.mismatches());

    ar = new AlignmentResult(DnaUtils.encodeString("acgtacgta"), ActionsHelper.build("X=======X", 0, 0), DnaUtils.encodeString("atgttcgta"));
    assertEquals(2, ar.mismatches());

    ar = new AlignmentResult(DnaUtils.encodeString("acgtacgta"), ActionsHelper.build("=X==X====", 0, 0), DnaUtils.encodeString("atgttcgta"));
    assertEquals(2, ar.mismatches());

    ar = new AlignmentResult(DnaUtils.encodeString("acgcgta"), ActionsHelper.build("===DD====", 0, 3), DnaUtils.encodeString("acgtacgta"));
    assertEquals(2, ar.mismatches());

    ar = new AlignmentResult(DnaUtils.encodeString("acgcgta"), ActionsHelper.build("===DD====", 0, 3), DnaUtils.encodeString("acgt"));
    assertEquals(1, ar.mismatches());

    ar = new AlignmentResult(DnaUtils.encodeString("acgcgta"), ActionsHelper.build("===DD====", 0, 3), DnaUtils.encodeString("acgtacgt"));
    assertEquals(2, ar.mismatches());

    ar = new AlignmentResult(DnaUtils.encodeString("acgtacgta"), ActionsHelper.build("===II====", 0, 3), DnaUtils.encodeString("acgcgta"));
    assertEquals(2, ar.mismatches());

    ar = new AlignmentResult(DnaUtils.encodeString("acgtacgta"), ActionsHelper.build("===II====", 0, 3), DnaUtils.encodeString("acg"));
    assertEquals(0, ar.mismatches());

    ar = new AlignmentResult(DnaUtils.encodeString("acgtacgta"), ActionsHelper.build("===II====", 0, 3), DnaUtils.encodeString("acgcgt"));
    assertEquals(2, ar.mismatches());
  }

  public void testMismatchesSoftClip() {
    AlignmentResult ar = new AlignmentResult(DnaUtils.encodeString("acgtacgta"), ActionsHelper.build("=========", 0, 0), DnaUtils.encodeString("acg"));
    assertEquals(0, ar.mismatches());
    ar = new AlignmentResult(DnaUtils.encodeString("acgtacgta"), ActionsHelper.build("SSSSSS=X=", 0, 0), DnaUtils.encodeString("gca"));
    assertEquals(1, ar.mismatches());
    ar = new AlignmentResult(DnaUtils.encodeString("acgtacgta"), ActionsHelper.build("SSSSSS=X=", 0, 0), DnaUtils.encodeString("gcatttttt"));
    assertEquals(1, ar.mismatches());
    ar = new AlignmentResult(DnaUtils.encodeString("acgtacgta"), ActionsHelper.build("=X=SSSSSS", 0, 0), DnaUtils.encodeString("aag"));
    assertEquals(1, ar.mismatches());
  }

  public void testReadString() {
    AlignmentResult ar = new AlignmentResult(DnaUtils.encodeString("acgtacgta"), ActionsHelper.build("=========", 0, 0), "".getBytes());
    ar.setIdentifyingInfo(true, false);
    assertEquals("ACGTACGTA", ar.readString());

    ar = new AlignmentResult(DnaUtils.encodeString("acgtacg"), ActionsHelper.build("==BB==R==", 0, 3), "".getBytes());
    ar.setIdentifyingInfo(true, false);
    assertEquals("ACACG", ar.readString());
    ar = new AlignmentResult(DnaUtils.encodeString("acgtacgt"), ActionsHelper.build("==BB=I==T=", 0, 3), "".getBytes());
    ar.setIdentifyingInfo(true, false);
    assertEquals("ACCGT", ar.readString());

    ar = new AlignmentResult(DnaUtils.encodeString("acgtacg"), ActionsHelper.build("===BB====", 0, 3), "".getBytes());
    ar.setIdentifyingInfo(false, false);
    assertEquals("ATACG", ar.readString());

    ar = new AlignmentResult(DnaUtils.encodeString("acgtacgt"), ActionsHelper.build("==I=BBN==N==", 0, 3), "".getBytes());
    ar.setIdentifyingInfo(false, false);
    assertEquals("AACGT", ar.readString());

    ar = new AlignmentResult(DnaUtils.encodeString("acgtacg"), ActionsHelper.build("==D=BB====", 0, 3), "".getBytes());
    ar.setIdentifyingInfo(false, false);
    assertEquals("ACTACG", ar.readString());

    ar = new AlignmentResult(DnaUtils.encodeString("acgtacg"), ActionsHelper.build("==N=BB====", 0, 3), "".getBytes());
    ar.setIdentifyingInfo(false, false);
    assertEquals("ACTACG", ar.readString());

    final int[] actions = ActionsHelper.build("==I=========D==", 0, 3);
    ActionsHelper.softClip(actions, true, 3, 0);
    ActionsHelper.softClip(actions, false, 2, 1);
    ar = new AlignmentResult(DnaUtils.encodeString("AATAAAAAAAAAGG"), actions, "".getBytes());
    ar.setIdentifyingInfo(false, false);
    assertEquals("AATAAAAAAAAAGG", ar.readString());
  }

  public void testCGAttrNew() {
    final byte[] qualityScores = new byte[35];
    for (int i = 0; i < qualityScores.length; ++i) {
      qualityScores[i] = (byte) i;
    }
    String read, acts, tmpl;
    AlignmentResult ar;
    BinaryTempFileRecord sr;

    //First arm, !rc
    read = "tcatg   a ggcacatcg gcattacagc       cagactgtcg";
    acts = "=====BBB=D=X======= ========== NNNN  ==========".replaceAll(" ", "");
    tmpl = "tcatg      acacatcg gcattacagc aggc  cagactgtcg";
    ar = new AlignmentResult(getSequence(read), ActionsHelper.build(acts, 0, 0), getSequence(tmpl));
    ar.setIdentifyingInfo(true, false);
    ar.setRemainingOutput(1, 0);
//    assertEquals("tcatggcacatcggcattacagc....cagactgtcg\ttcatgacacatcggcattacagcaggccagactgtcg\t||||| |||||||||||||||||    ||||||||||", ar.tabularString());

    sr = ar.toRecord(true, null, 0, false, false);
    assertEquals("5=1X17=4N10=", new String(sr.getCigarString()));
    assertEquals("5=3B1=1D1=1X17=4N10=", new String(sr.getSuperCigarString()));
    assertEquals("G", new String(sr.getReadDeltaString()));
    assertEquals("TCATGGCACATCGGCATTACAGCCAGACTGTCG", new String(sr.getCgReadString()));
    assertEquals(33, sr.getCgReadString().length);
    // first arm, !rc again, but with no errors.
    read = "tcatg   atgcacatcg gcattacagc       cagactgtcg";
    acts = "=====BBB========== ========== NNNN  ==========".replaceAll(" ", "");
    tmpl = "tcatg      cacatcg gcattacagc aggc  cagactgtcg";
    ar = new AlignmentResult(getSequence(read), ActionsHelper.build(acts, 0, 0), getSequence(tmpl));
    ar.setIdentifyingInfo(true, false);
    ar.setRemainingOutput(1, 0);
    sr = ar.toRecord(true, null, 0, false, false);
    assertEquals("22=4N10=", new String(sr.getCigarString()));
    assertEquals("5=3B20=4N10=", new String(sr.getSuperCigarString()));
    assertEquals(0, sr.getReadDeltaString().length);
    assertEquals("TCATGCACATCGGCATTACAGCCAGACTGTCG", new String(sr.getCgReadString()));

    //First arm, rc
    read = "tattt   ttcgaggact tattttaatt       ctcaaacgcg";
    //tmprc tattt        ggact tattttaatt aattg ctcaaacgct
    acts = "=====BBB==IIX===== ========== NNNNN =========X".replaceAll(" ", "");
    tmpl = "agcgtttgag caatt aattaaaata agtcc aaata         agagagagagagaga"; // is reverse complement compared to the read
    ar = new AlignmentResult(getSequence(read), ActionsHelper.build(acts, 0, 0), getSequence(tmpl));
    //====D===================================          ====D===================================
    ar.setIdentifyingInfo(true, true);
    ar.setRemainingOutput(1, 0);
//    assertEquals("tatttggacttattttaatt.....ctcaaacgcg\ttatttggacttattttaattaattgctcaaacgct\t||||||||||||||||||||     ||||||||| ", ar.tabularString());

    assertEquals("TATTT" + "GGACT" + "TATTTTAATT" + "CTCAAACGCG", ar.readString());
    sr = ar.toRecord(true, null, 0, false, false);
    assertEquals("1X9=5N20=", new String(sr.getCigarString()));
    assertEquals("1X9=5N15=1X2I2=3B5=", new String(sr.getSuperCigarString()));
    assertEquals("CTCG", new String(sr.getReadDeltaString()));
    assertEquals("CGCGTTTGAGAATTAAAATAAGTCCAAATA", new String(sr.getCgReadString()));
    assertEquals(30, sr.getCgReadString().length);


    //second arm, rc
    read = "ctgctgaccc          gacctactgg cggactcgag  agtta";
    acts = "========== NNNNNNNN ========== ========X=BBX====".replaceAll(" ", "");
    tmpl = "ctgctgaccc agagccgg gacctactgg cggactcg    ggtta";
    ar = new AlignmentResult(getSequence(read), ActionsHelper.build(acts, 0, 0), getSequence(tmpl));
    ar.setIdentifyingInfo(false, true);
    ar.setRemainingOutput(1, 0);

    sr = ar.toRecord(true, null, 0, false, false);
    assertEquals("4=1X18=8N10=", new String(sr.getCigarString()));
    assertEquals("4=1X2B1=1X18=8N10=", new String(sr.getSuperCigarString()));
    assertEquals("TT", new String(sr.getReadDeltaString()));
    assertEquals("TAACTCGAGTCCGCCAGTAGGTCGGGTCAGCAG", new String(sr.getCgReadString()));
    assertEquals(33, sr.getCgReadString().length);

    //second arm, !rc
    read = "ttccggtctt         ctctctatcc cgctggtttt  acagg";
    acts = "========== NNNNNNN ========== ========XXBB=I===".replaceAll(" ", "");
    tmpl = "ttccggtctt tgtcttt ctctctatcc cgctggtt    a agg                  attatatatat";
    ar = new AlignmentResult(getSequence(read), ActionsHelper.build(acts, 0, 0), getSequence(tmpl));
    ar.setIdentifyingInfo(false, false);
    ar.setRemainingOutput(1, 0);

    sr = ar.toRecord(true, null, 0, false, false);
    assertEquals("10=7N19=1I3=", new String(sr.getCigarString()));
    assertEquals("10=7N18=2X2B1=1I3=", new String(sr.getSuperCigarString()));
    assertEquals("TTC", new String(sr.getReadDeltaString()));
    assertEquals("TTCCGGTCTTCTCTCTATCCCGCTGGTTACAGG", new String(sr.getCgReadString()));
    assertEquals(33, sr.getCgReadString().length);

    //second arm, !rc  (again, but with the whole overlap counteracted by insertions)
    read = "ttccggtctt         ctctctatcc cgctggtttt    atagg";
    acts = "========== NNNNNNN ========== ==========DDBB=I===".replaceAll(" ", "");
    tmpl = "ttccggtctt tgtcttt ctctctatcc cgctggtttt    a agg                              atatatatat";
    ar = new AlignmentResult(getSequence(read), ActionsHelper.build(acts, 0, 0), getSequence(tmpl));
    ar.setIdentifyingInfo(false, false);
    ar.setRemainingOutput(1, 0);

    sr = ar.toRecord(true, null, 5, false, false);
    assertEquals("10=7N21=1I3=", new String(sr.getCigarString()));
    assertNotNull(sr.getSuperCigarString());
    assertEquals("10=7N20=2D2B1=1I3=", new String(sr.getSuperCigarString()));
    assertNotNull(sr.getReadDeltaString());
    assertEquals("T", new String(sr.getReadDeltaString()));
    assertEquals("TTCCGGTCTTCTCTCTATCCCGCTGGTTTTATAGG", new String(sr.getCgReadString()));
    assertEquals(35, sr.getCgReadString().length);
    assertEquals(6, sr.getStartPosition());
  }

  //this is \n and not LS as stupid SAM uses "\n" need to fix this
  private static final String SAM_SEQ = "ggataggagagtacagataattagatgttttgccc";
  private static final byte[] SAM_QUALITY = new byte[SAM_SEQ.length()];

  static {
    for (int k = 0; k < SAM_QUALITY.length; ++k) {
      SAM_QUALITY[k] = (byte) k;
    }
  }

  private static String bcarToSamString(BinaryTempFileRecord rec) {
    //0#0#toxic#1#255#35=#*#0#0#*#*#AS:i:1#NM:i:0
    return rec.getReadId() + "\t"
        + (rec.getSamFlags() & 0xFF)
        + "\t" + rec.getReferenceId()
        + "\t" + rec.getStartPosition()
        + "\t*\t" + new String(rec.getCigarString())
        + "\t=\t" + rec.getMatePosition()
        + "\t" + rec.getTemplateLength()
        + "\t*\t*\tAS:i:"
        + rec.getAlignmentScore()
        + "\tNM:i:" + rec.getNumberMismatches();
  }

  private void runSAMPairedEnd(AlignmentResult mate, int pos, String complement, String expected, boolean first) throws Exception {
    final AlignmentResult ar = new AlignmentResult(getSequence(SAM_SEQ), ActionsHelper.build("==========NNNNN=========================", pos - 1, 1), getSequence(SAM_SEQ + "aaaaaaa"));
    ar.setIdentifyingInfo(first, complement.equals("R"));
    ar.mReferenceId = 7;
    final File f = File.createTempFile("junit", "sam");
    try {
      final BinaryTempFileRecord sam = ar.toRecord(mate != null, mate, 0, false, false);

      final String x = bcarToSamString(sam).replace('\t', '#');
      assertEquals(x, expected, x);

      assertEquals(complement.equals("R"), sam.isReverseStrand());
      assertEquals((sam.getSamFlags() & SamBamConstants.SAM_READ_IS_PAIRED) == 0, mate == null);

      if (mate != null) {
        assertTrue(sam.isReadPaired());
        assertEquals(mate.getStart(), sam.getMatePosition() - 1);
        assertEquals(first, sam.isFirstOfPair());
        assertEquals(!first, sam.isSecondOfPair());
      } else {
        assertFalse(sam.isReadPaired());
        assertEquals(0, sam.getMatePosition());
        assertFalse(sam.isFirstOfPair());
        assertFalse(sam.isSecondOfPair());
      }
    } finally {
      assertTrue(f.delete());
    }
  }

  private void runSAM(final String complement, final String expected) throws Exception {
    runSAMPairedEnd(new AlignmentResult(null, ActionsHelper.build("", 0, 0), null), 1, complement, expected, false);
  }

  public void testSamTrueSingleEnd() throws Exception {
    runSAMPairedEnd(null, 0, "R", "0#16#7#1#*#1S24=5N10=#=#0#0#*#*#AS:i:1#NM:i:0", false);
  }

  public void testSam() throws Exception {
    runSAM("F", "0#131#7#1#*#10=5N25=#=#1#-40#*#*#AS:i:1#NM:i:0");
  }

  public void testSAMRC() throws Exception {
    runSAM("R", "0#147#7#1#*#25=5N10=#=#1#-40#*#*#AS:i:1#NM:i:0");
  }

  public void testSAMQ() throws Exception {
    runSAM("F", "0#131#7#1#*#10=5N25=#=#1#-40#*#*#AS:i:1#NM:i:0");
  }

  public void testSAMQRC() throws Exception {
    runSAM("R", "0#147#7#1#*#25=5N10=#=#1#-40#*#*#AS:i:1#NM:i:0");
  }

  public void testSAMPE() throws Exception {
    final AlignmentResult mate = new AlignmentResult(getSequence("ACGT"), ActionsHelper.build("", 12, 0), getSequence(SAM_SEQ));
    mate.mReferenceId = 6;
    mate.setIdentifyingInfo(false, true);
    runSAMPairedEnd(mate, 1, "F", "0#99#7#1#*#10=5N25=#=#13#40#*#*#AS:i:1#NM:i:0", true);
  }

  public void testSAMPEInsert() throws Exception {
    final AlignmentResult mate = new AlignmentResult(getSequence("ACNGT"), ActionsHelper.build("", 12, 4), getSequence(SAM_SEQ));
    mate.mReferenceId = 8;
    mate.setIdentifyingInfo(false, true);
    runSAMPairedEnd(mate, 1, "F", "0#163#7#1#*#10=5N25=#=#13#40#*#*#AS:i:1#NM:i:0", false);
  }

  public void testSAMPEDeletion() throws Exception {
    final AlignmentResult mate = new AlignmentResult(getSequence("ACGT"), ActionsHelper.build("", 12, 2), getSequence("ACNT"));
    mate.mReferenceId = 7;
    mate.setIdentifyingInfo(false, true);
    runSAMPairedEnd(mate, 1, "F", "0#163#7#1#*#10=5N25=#=#13#40#*#*#AS:i:1#NM:i:0", false);
  }

  public void testSAMPESoftLeftClip() throws Exception {
    final AlignmentResult mate = new AlignmentResult(getSequence("ACGT"), ActionsHelper.build("", 12, 2), getSequence(SAM_SEQ));
    mate.mReferenceId = 0;
    mate.setIdentifyingInfo(false, true);
    runSAMPairedEnd(mate, -2, "F", "0#163#7#1#*#3S7=5N25=#=#13#40#*#*#AS:i:1#NM:i:0", false);
  }

  public void testBug867() {
    final AlignmentResult ar = new AlignmentResult(getSequence("taaataatggcaatatctgcaggga     aaactataag"), ActionsHelper.build("", 0, 0), getSequence("taaNNNatggcaatatctgcagggaacaggaaactataag"));
    assertEquals(0, ar.mismatches());
  }

  public void testCgTripleInsert() {
    final String read =    "tagacaaatg        ttacaagaccacaggaggggaa".replaceAll(" ", "");
    final String temp =    "tagacaaatgtgactggattacaagccacaggagggggaaa";
    final String actions = "==========NNNNNNND=======I============X=";
    final AlignmentResult alignment2 = new AlignmentResult(DnaUtils.encodeString(read), ActionsHelper.build(actions, 0, 0), DnaUtils.encodeString(temp));

    alignment2.setIdentifyingInfo(false, false);
    alignment2.setRemainingOutput(-2, 3);

    assertEquals(0, alignment2.getScore());
    assertFalse(alignment2.isFirst());
    assertEquals(0, alignment2.getStart());
    assertEquals(-2, alignment2.getReadId());
    assertEquals(3, alignment2.getReferenceId());
    assertEquals(actions, alignment2.getActionsString());
    assertEquals(30, alignment2.getMatchCount());
    assertEquals(1, alignment2.getDeletionsFromReadCount());
    assertFalse(alignment2.isReverse());
    assertEquals("TAGACAAATGTTACAAGACCACAGGAGGGGAA", alignment2.readString());

    final String cigar2 = alignment2.getCigarString(false, false);
    assertEquals("10=7N1D7=1I12=1X1=", cigar2);

  }
  public void testCgTripleInsertRev() {
    final String read = "                           tagacaaatg        ttacaagaccacaggaggggaa".replaceAll(" ", "");
    final String temp = DnaUtils.reverseComplement("tagacaaatgtgactggattacaag ccacaggagggggaaa".replaceAll(" ", ""));

    final String actions = "==========NNNNNNND=======I============X=";
    final AlignmentResult alignment2 = new AlignmentResult(DnaUtils.encodeString(read), ActionsHelper.build(actions, 2, 0), DnaUtils.encodeString(temp));

    alignment2.setIdentifyingInfo(false, true);
    alignment2.setRemainingOutput(-2, 1);

    assertEquals(0, alignment2.getScore());
    assertFalse(alignment2.isFirst());
    assertEquals(2, alignment2.getStart());
    assertEquals(-2, alignment2.getReadId());
    assertEquals(1, alignment2.getReferenceId());
    assertEquals(actions, alignment2.getActionsString());
    assertEquals(30, alignment2.getMatchCount());
    assertEquals(1, alignment2.getDeletionsFromReadCount());
    assertTrue(alignment2.isReverse());
    assertEquals("TAGACAAATGTTACAAGACCACAGGAGGGGAA", alignment2.readString());

    final String cigar2 = alignment2.getCigarString(true, false);
    assertEquals("1=1X12=1I7=1D7N10=", cigar2);
  }

  public void testCompare() {
    final AlignmentResult ar1 = new AlignmentResult(new byte[0], new int[3], new byte[0]);
    ar1.setIdentifyingInfo(false, true);
    ar1.setRemainingOutput(-2, 0);

    final AlignmentResult ar2 = new AlignmentResult(new byte[0], new int[3], new byte[0]);
    ar2.setIdentifyingInfo(false, true);
    ar2.setRemainingOutput(-2, 0);

    assertTrue(ar1.equals(ar2));

    ar2.setIdentifyingInfo(false, false);
    assertFalse(ar1.equals(ar2));

    ar2.setIdentifyingInfo(false, true);
    ar2.setRemainingOutput(-2, 1);
    assertFalse(ar1.equals(ar2));

    final AlignmentResult ar3 = new AlignmentResult(new byte[0], ActionsHelper.build("", 2, 3), new byte[0]);
    ar3.setIdentifyingInfo(false, true);
    ar3.setRemainingOutput(-2, 0);

    assertFalse(ar1.equals(ar3));
    assertNotNull(ar1);
  }

  public void testOverlapCigarConsistency() {
    final byte[] read = DnaUtils.encodeString("GTACTGTGTGCGTCCTTGACATCTA     ACCGGCGCCT".replaceAll(" ", ""));
    final byte[] tmpl = DnaUtils.encodeString("GTAGT  GTGCGTCCTTGACATCTAAGGATATCGGCGCCTGA".replaceAll(" ", ""));
    final int[] actions = ActionsHelper.build("===X=BB====================NNNNN=X========", 0, 3);
    final AlignmentResult ar = new AlignmentResult(read, actions, tmpl);
    ar.setIdentifyingInfo(true, false);
    //System.out.println("cigar:" + ar.getCigarString(1, false) + " readString:" + ar.readString());
    final SAMRecord samrec = new SAMRecord(null);
    samrec.setCigarString(ar.getCigarString(false, false));
    samrec.setReadString(ar.readString());
    samrec.setAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES, 2);
    samrec.setReadNegativeStrandFlag(false);
    samrec.setFirstOfPairFlag(true);
    samrec.setAlignmentStart(1);
    samrec.setFlags(67);
    final NgsParams params = new NgsParamsBuilder().gapOpenPenalty(EditDistanceFactory.DEFAULT_GAP_OPEN_PENALTY)
        .gapExtendPenalty(EditDistanceFactory.DEFAULT_GAP_EXTEND_PENALTY)
        .substitutionPenalty(EditDistanceFactory.DEFAULT_SUBSTITUTION_PENALTY)
        .unknownsPenalty(0).create();
    try (MemoryPrintStream mps = new MemoryPrintStream()) {
      final SamValidator sv = new SamValidator(mps.printStream(), mps.printStream(), true, false, false, false, params, true);
      assertEquals(2 * EditDistanceFactory.DEFAULT_SUBSTITUTION_PENALTY, sv.isAtExpectedRef(tmpl, samrec, null));
    }
  }
  public void testOverlapCigarConsistencyRev() {
    final byte[] read = DnaUtils.encodeString(
                   DnaUtils.reverseComplement("TAGGCGGGTTGCCAA TTAACTTGTA      GTCCTTGACA".replaceAll(" ", "")));
    final byte[] tmpl = DnaUtils.encodeString("TAGGG   TGGCCAA TTAACTTGTAGTGTGCGTCCTTGACA".replaceAll(" ", ""));
    final int[] actions = ActionsHelper.build("==========NNNNNN===============X====BBBX====", 0, 3);
    final AlignmentResult ar = new AlignmentResult(read, actions, tmpl);
    ar.setIdentifyingInfo(false, true);
    //System.out.println("cigar:" + ar.getCigarString(1, true) + " readString:" + ar.readString());
    final SAMRecord samrec = new SAMRecord(null);
    samrec.setCigarString(ar.getCigarString(true, false));
    samrec.setReadString(DnaUtils.reverseComplement(ar.readString()));
    samrec.setAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES, 2);
    samrec.setAlignmentStart(1);
    samrec.setFlags(179);
    final NgsParams params = new NgsParamsBuilder().gapOpenPenalty(EditDistanceFactory.DEFAULT_GAP_OPEN_PENALTY)
        .gapExtendPenalty(EditDistanceFactory.DEFAULT_GAP_EXTEND_PENALTY)
        .substitutionPenalty(EditDistanceFactory.DEFAULT_SUBSTITUTION_PENALTY)
        .unknownsPenalty(0).create();
    try (MemoryPrintStream mps = new MemoryPrintStream()) {
      final SamValidator sv = new SamValidator(mps.printStream(), mps.printStream(), true, false, false, false, params, false);
      assertEquals(2 * EditDistanceFactory.DEFAULT_SUBSTITUTION_PENALTY, sv.isAtExpectedRef(tmpl, samrec, null));
    }
  }

  public void testInvalidSoftClipping() {
    final byte[] read = DnaUtils.encodeString("GAGGGTTAGG     GTGAGGGTTTGGGTTAGGGTATTAG".replaceAll(" ", ""));
    byte[] tmpl = DnaUtils.encodeString("GAGGGTTAGGGTTAGGGTGAGGGTTAGGGTTAGGG".replaceAll(" ", ""));
    int[] actions = ActionsHelper.build("==========NNNNNN=========X==========B=====", 0, 3);
    AlignmentResult ar = new AlignmentResult(read, actions, tmpl);
    assertEquals("10=6N9=1X9=5S", ar.getCigarString(false, false));

              //                            GAGGGTTAGG     .GTGAGGGTTTGGGTTAGGGTATTAG
    tmpl = DnaUtils.encodeString("AAAAAAAAAAGAGGGTTAGGGTTAGGGTGAGGGTTAGGGTTAGGG".replaceAll(" ", ""));
    actions = ActionsHelper.build("==========NNNNNN=========X==========B=====", 10, 3);
    ar = new AlignmentResult(read, actions, tmpl);
    assertEquals("10=6N9=1X9=5S", ar.getCigarString(false, false));
//    assertEquals("gagggttagg......gtgagggtttgggttagggattag\tgagggttagggttagggtgagggttagggttagggnnnnn\t||||||||||      ||||||||| |||||||||     ", ar.tabularString());
  }
}
