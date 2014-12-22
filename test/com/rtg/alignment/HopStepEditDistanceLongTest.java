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

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;

import com.rtg.mode.DnaUtils;
import com.rtg.ngs.NgsParams;
import com.rtg.ngs.NgsParamsBuilder;
import com.rtg.util.MaxShiftUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;


/**
 * Tests corresponding class
 */
public class HopStepEditDistanceLongTest extends AbstractUnidirectionalEditDistanceTest {

  private HopStepEditDistanceLong mHopStep;
  private int mMaxScore;

  @Override
  protected UnidirectionalEditDistance getEditDistanceInstance(int gapOpenPenalty, int gapExtendPenalty, int subsPenalty, int unknownsPenalty) {
    final NgsParams params = new NgsParamsBuilder().gapOpenPenalty(gapOpenPenalty).gapExtendPenalty(gapExtendPenalty).substitutionPenalty(subsPenalty).unknownsPenalty(unknownsPenalty).create();
    final HopStepEditDistanceLong ed = new HopStepEditDistanceLong(params);
    ed.setReadLengthCutoff(10);
    return ed;
  }

  @Override
  public void setUp() throws Exception {
    super.setUp();
    mMaxScore = 50;
    final NgsParams params = new NgsParamsBuilder().gapOpenPenalty(1).gapExtendPenalty(1).substitutionPenalty(1).unknownsPenalty(0).create();
    mHopStep = new HopStepEditDistanceLong(params);
    mHopStep.setReadLengthCutoff(10);
  }

  @Override
  public void tearDown() throws Exception {
    super.tearDown();
    mHopStep = null;
  }

  public int[] align(String read, String template, int zeroBasedPos, int mismatches, boolean resnull, int maxShift) {
    final byte[] encodedRead = DnaUtils.encodeString(read.replaceAll(" ", ""));
    final byte[] encodedTemplate = DnaUtils.encodeString(template.replaceAll(" ", ""));
    final int[] actions = mHopStep.calculateEditDistance(encodedRead, encodedRead.length, encodedTemplate, zeroBasedPos, mMaxScore, maxShift, true);
    mHopStep.logStats();
    mHopStep.globalIntegrity();
    if (resnull) {
      if (actions != null) {
        System.err.println("as=" + actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
      }
      assertNull(actions);
    } else {
      assertNotNull(actions);
      mHopStep.validateResult(zeroBasedPos, Integer.MAX_VALUE);
      //System.out.println(ActionsHelper.toString(actions));
      final AlignmentResult result = new AlignmentResult(encodedRead, actions, encodedTemplate);
      assertEquals(mismatches, result.mismatches());
    }
    return actions;
  }

  public void testalign1() {
    align("aTgtacgtacgtacgtacgtacgtacgaacgtacgt",
          "aCgtacgtacgtacgtacgtacgtacgaacgtacgt", 0, 1, false, 4);
  }

  // 2 nearby SNPs should be aligned using full aligner
  public void testalign2() {
    align("atatatatatataaatCtCtatatatatataaaa",
          "atatatatatataaatatatatatatatataaaa", 0, 2, false, 4);
  }

  // lots of Ns, even in start region.
  public void testalign3() {
    align("atatatatatatatatatatatatatatataaaa",
          "NtatatatatatNNNNNNNNatatatNtataNaa", 0, 0, true, 4);
  }

  public void testalign3Mismatch2() {
    align("atatatatatatatatatatatatatatataaaa",
          "CtatatatatatatatatatatatatCtataCaa", 0, 3, true, 4);
  }//was:  XXXXX=========================XX=X

  public void testalign4() {
    align("aaagaaaaaaaaa",
          "aaagaaaaaaaaa", 0, 0, false, 4);
  }

  // fails to initially synchronize
  public void testalign5() {
    align("atatatatatatatatatatatatatatatatatatatatat",
          "atatatatatatatatatatatatatatatatatatatatat", 0, 0, true, 4);
  }

  // has Ns before SNP
  public void testalign6() {
    align("atatatatatatatatatataNNtttatAtatatatataaat",
          "atataNatatatatatatatatatttatGtatatatataNat", 0, 1, true, 4);
  }

  // too short to synchronize
  public void testalign7() {
    align("ata atat",
          "atatatat", -1, 1, true, 4);
  }

  // cannot resynchronize after indel, due to a repeat region
  public void testalign8() {
    align("atatatatatatatatAtatcgtttta",
          "atatatatatatatat tatcgtttta", 0, 1, true, 4);
  }

  public void testalign9() {
    align("actgactagctagctagctgactgaCtgactgactga",
         "aactgactagctagctagctgactga tgactgactga", 0, 1, false, 4);
  }

  // has to synchronize quite a long way down from the end
  public void testalign19() {
    align("actgactagctagctagctgacctgactgactgactga",
         "tactgactagctagctagctgacctgactgactgactga", 1, 0, false, 4);
  }

  // "CC" deleted in middle of read
  public void testalign20() {
    align("actgactagctagctagctgac  tgactgcctgactga",
         "tactgactagctagctagctgacCCtgactgcctgactga", 0, 2, false, 4);
  }

  public void testalign20b() {
    align(
                            "actgactagctagctagctga  ctgactgactgactgg",
       "nnnnnnnnnnnnnnnnnnnntactgactagctagctagctgaCCctgactgactgactgg", 22, 2, false, 4);
  }

  public void testalign20OffTemplate() {
    align(
                            "actgactagctagctagctga  ctgattgactgactga",
       "nnnnnnnnnnnnnnnnnnnntactgactagctagctagctgaCCctgattgact", 22, 2, true, 4);
  }

  public void testalign20OffTemplateLeft() {
    align(
       "actgactagctagctagctga  ctgattgactgactga",
       " ctgactagctagctagctgaCCctgattgactgactg", 0, 2, true, 4);
  }

  public void testalign21() {
    align("aCtgccccccccc",
          "agtgccccccccc", 0, 1, false, 4);
  }

  public void testalign22() {
    align("ACtgacacgacac",
          "ggtgacacgacac", 0, 2, false, 4);
  }

  public void testalign23() {
    align("AttttatttttttCtgccccccccc",
          "GttttatttttttGtgccccccccc", 0, 2, false, 4);
  }

  // a SNP near the end.
  public void testalign24() {
    align("AtttttcgcgttttCtg",
          "gtttttcgcgttttgtg", 0, 2, false, 4);
  }

  public void testalign25() {
    final int[] actions = align("AtttttgttttCtgccccccccc",
          "gtttttgttttgtgccccccccc", 0, 2, false, 4);
    assertNotNull(actions);
  }

  public void testalign25offset10() {
    final int[] actions =
    align(
        "AAttttgttttCtgccccccccc",
        "gtttttgttttgtgcccccccccaaaaaaaaaa", 4, 3, false, 4);
    assertEquals(3, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
    assertEquals(0, actions[ActionsHelper.TEMPLATE_START_INDEX]);
  }

  public void testalign25offsetneg10() {
    align("          Atttttgttttctgccccccccc",
          "ccccccccccgtttttgttttgtgccccccccc", 6, 2, false, 4);
  }

  public void testalign26() {
    final int[] actions =
      align("ctgggttgagctTgaagtggtcgtgacaggcgcacggtttcccaacgacggaaggttaaAgctaccattggtagtgaagcaggcgtattgactcggctgg",
           "cctgggttgagct gaagtggtcgtgacaggcgcacggtttcccaacgacggaaggttaaTgctaccattggtagtgaagcaggcgtattgactcggctgg",
           0, 2, false, 4);
    assertEquals(3, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
    assertEquals(1, actions[ActionsHelper.TEMPLATE_START_INDEX]);
  }

  public void testalign27() {
    final int[] actions =
    align("ctcccctgtagccttgctagcggttcgaaggaCtgtgcttttgcTGGttccttatccgtacagttggtGtggcgcagcttcttacgttagactgaccatC",
        "tactcccctgtagccttgctagcggttcgaaggaGtgtgcttttgc   ttccttatccgtacagttggtCtggcgcagcttcttacgttagactgaccatT", 0, 6, false, 4);
    final String actString =
          "================================X===========III=====================X==============================X";
    assertEquals(actString, ActionsHelper.toString(actions));
    assertEquals(7, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
  }

  public void testalign28() {
    final int[] actions =
    align("ctcccctgtagccttgctagcggttcgaaggaCtgtgcttttgc   gttccttatccgtacagttggtGtggcgcagcttcttacgttagactgCccatC",
          "ctcccctgtagccttgctagcggttcgaaggaGtgtgcttttgcTGAgttccttatccgtacagttggtCtggcgcagcttcttacgttagactgAccatT", 0, 7, false, 4);
    final String actString =
          "================================X===========DDD======================X=========================X====X";
    assertEquals(actString, ActionsHelper.toString(actions));
    assertEquals(8, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
  }

  // MNP of length 3 at position 7
  public void testalign29() {
    align("ATATGAAGACCAGCAATCTATCTTGATCCACAGGCGCGTCGCTCAAATACTCGGCTAGATGCGCAGTCACGCCCTATATGAAGAAAGATGCTCCTACGAA",
          "ATATGAAagaCAGCAATCTATCTTGATCCACAGGCGCGTCGCTCAAATACTCGGCTAGATGCGTAGTCACGCCCTATATGAAGAAAGAcGCTCCTACGAA",
          0, 5, false, 4);
  }

  public void testalign30() {
    align("CCGTTCCGCTGCAGGGCCCTGGCATGTATGCGTTGTACATCACATAAGTTGGTCTAGTCTTAGGCCAATTGCCAATCCGGTGGTCTCGACATACCGCGAA",
          "atatatatatatatatatatatatatatatCAGTTCCGCTGCAGGGCCCTGGCATGTATGCGTTGTACCTCACATAAGTTGGTCTAGTCTTAGGCCAATTGCCAATCCGGTGGTCTCGACATACCGCGAA",
          30, 2, false, 4);
  }

  // two many indels too close
  public void testalign31() {
    align("agcatcgtagcatgtagcTatatttggcattcAgcatggctacacattctaTACgctagccccttttgctag  ctgatcgatcgatcgatcgatcc",
        "aagcatcgtagcatgtagc atatttggcattc gcatggctacacattcta CAgctagccccttttgctagCCctgatcgatcgatcgatcgatcc",
        1, 6, true, 4);
  }
  public void testalign31b() {
    align("agcatcgtagcatgtagcTatatgcattcAgcatggcctacacattctaTACgctagccccttttgctag  ctgatcgatcgatcgatcgatcc",
         "aagcatcgtagcatgtagc atatgcattc gcatggcctacacattcta CAgctagccccttttgctagCCctgatcgatcgatcgatcgatcc",
         1, 6, false, 4);
  }

  // test that it can move the start position forwards/backwards by MAX_INDEL + 1 positions with read length above 100.
  public void testalignkurt1() {
    final String read =           "cacggatcagctctacgttccttACGTTgccacggatcagctctacgttccttgaagcgcGtccgccgaacggcagcaatcactacttacgacatgtctaccc";
    final String tmpl = "ttcaagtgtccacggatcagctctacgttccttgaagcgccacggatcagctctacgttccttgaagcgcatccgccgaacggcagcaatcactacttacgacatgtctaccctcgactcgac";
    int[] actions = align(read, tmpl, 4, 6, true, 5);
    actions = align(read, tmpl, 5, 6, false, 5);
    assertEquals(6, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
    assertEquals(10, actions[ActionsHelper.TEMPLATE_START_INDEX]);
    actions = align(read, tmpl, 15, 6, false, 5);
    assertEquals(6, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
    assertEquals(10, actions[ActionsHelper.TEMPLATE_START_INDEX]);
    actions = align(read, tmpl, 16, 6, true, 5);
    assertNull(actions);
  }

  // the beginning of the read is too nasty
  public void testalignStartTooHard() {
    mMaxScore = 5;
    align("tgggaatgcttgacaaaaacaagcaatggggaaaggattccctatttaataaatg",
          "ATCTTTGAGCCgaGaaaaacaagcaatggggaaaggattccctatttaataaatg",
          0, 0, true, 4);
  }

  public void testalignMiddleTooHard() {
    mMaxScore = 5;
    align("tgggaatgcctgaCaaaaacaagcaatACCCCCCTCaaaggattccctatttaataaatg",
          "tgggaatgcctgaGaaaaacaagcaat ggggggggaaaggattccctatttaataaatg",
          0, 0, true, 4);
  }

  public void testalignEndTooHard() {
    mMaxScore = 5;
    align("tgggaatgcctgaCaaaaacaagcaataggggaaaggattccctatttaGaCCCCCC",
          "tgggaatgcctgaGaaaaacaagcaataggggaaaggattccctattta ataattg",
          0, 0, true, 4);
  }

  public void testalignMiddleShiftsTooFar() {
    align("tgggaatgcctgac  aaaaacaagcaat   aaacaggattccGGctca tttaataattg",
          "tgggaatgcctgacTTaaaaacaagcaatGGGaaacaggattcc  ctcaGtttaataattg",
          0, 6, true, 4);
  }

  /** This one was stopping too early, because 2 inserts are close together (around position 710) */
  public void test1000() {
    final int[] actions =
      align("actgcgaagctccacctgcacGacaggcctgtttcttgtctctttcttcttttctttcctttccttttaaaaatactttctagatttaatcatttttgagtgtgcagtCggaacatttaggaagatgccagacacaccagcaccgcccaggaacctgctgtcttggccggacacaccagcacagccaggaaactgtggtcacgacccgcggggcagaagcgagaaccttggggccggccccacccctggccgtgcctccctggcctggccctggggcctccttgtctttcctcccggttcaccagccaggcacacgcatccccaggcacacgcgtccccaggcacacgcgtccccTaggcacacgcgtccccaggctgcactggcttcgagctgtgtgttcctggagttgttgtctcacctggtttCcatagcagatttctctgatttctttccctgcacgtgtttctttaaggcccatcctgtaaccgtcggtagctggggctccccatgtccatggctgcacGgggacccgcttctggcgagaccacggggctttatccctttatccctcctgCTccgactgtggcacgttggtctccgagggctcaggtggtctgcaggctcgtcccacccatgcggggctcagctctaCgcctgctaactcctcacaagcccatccctggcaggggccccgcgagtccctctgcgcccgtccctTctgTcgcccgtctgtgagcGttgctggctgcgccgtcctctcttgcctctgtgaagcgtgtgcaggtctcttccgctctcttctgggccgcccactcttcctcactgacgggcaggagcctgtgtttgttcctcacatcagcccgagcggcagagacttt cctttcgctgagtcagctttggAaacccaaagAAtgtgatactggggtgcctcgtgcttctagggatggtttaaggtacccgtcctcatgccacatccacaccatgtctcctccacgaggc",
  "aaaaacctctactgcgaagctccacctgcacCacaggcctgtttcttgtctctttcttcttttctttcctttccttttaaaaatactttctagatttaatcatttttgagtgtgcagt ggaacatttaggaagatgccagacacaccagcaccgcccaggaacctgctgtcttggccggacacaccagcacagccaggaaactgtggtcacgacccgcggggcagaagcgagaaccttggggccggccccacccctggccgtgcctccctggcctggccctggggcctccttgtctttcctcccggttcaccagccaggcacacgcatccccaggcacacgcgtccccaggcacacgcgtcccc aggcacacgcgtccccaggctgcactggcttcgagctgtgtgttcctggagttgttgtctcacctggttt catagcagatttctctgatttctttccctgcacgtgtttctttaaggcccatcctgtaaccgtcggtagctggggctccccatgtccatggctgcac gggacccgcttctggcgagaccacggggctttatccctttatccctcctg  ccgactgtggcacgttggtctccgagggctcaggtggtctgcaggctcgtcccacccatgcggggctcagctctaagcctgctaactcctcacaagcccatccctggcaggggccccgcgagtccctctgcgcccgtccct ctg cgcccgtctgtgagcattgctggctgcgccgtcctctcttgcctctgtgaagcgtgtgcaggtctcttccgctctcttctgggccgcccactcttcctcactgacgggcaggagcctgtgtttgttcctcacatcagcccgagcggcagagactttCcctttcgctgagtcagctttgg aacccaaagTCtgtgatactggggtgcctcgtgcttctagggatggtttaaggtacccgtcctcatgccacatccacaccatgtctcctccacgaggctccccagtcaccccctgtgagca",
        5, 15, false, 10);
    final String aStr =
            "=====================X======================================================================================I======================================================================================================================================================================================================================================================I======================================================================I=================================================================================================I==================================================II===========================================================================X=================================================================I===I===============X============================================================================================================================================D======================I=========XX=======================================================================================";
    assertEquals(aStr, ActionsHelper.toString(actions));
    assertEquals(24, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
    assertEquals(10, actions[ActionsHelper.TEMPLATE_START_INDEX]);
  }

  /** Did have an incorrect alignment score. */
  public void test1000Score30() {
    final int[] actions =
      align(
                                  "aaggtattacaattgactccagaaaacaaaagggttgaagcacagggtggagggggctgaagccatgtggctcaatggccatgtg aagaaacaggccaacccttagtcccagttcactcagtgggggcctgatacagacggccccacagtgaaagacgggttaacaaaagtcaagagctcaggataagattcttctctggcgtaccaagatttttgttcttggtgccacacctcacccgatactgacaaaagcaaagcaggcaatccgagagaagtgggaaggctcaaacctatgtattaagaagatctggaaagattgggaattttaggcatgaagaagaaagagtcagaaaggacaagtcccacataaactgccatgtgggaggacatggaacttgtattctatggccactgagattagg tgtaggactaaaggagagaagctatagcaaagggattagggcttaatgagggaggactttcttctgcatgtactgcagcaagtaggagggagggggatagatgaccacttttcagggagttggaagagattcattcttccaagggatggtcgggcaagacagcctcaaggtccttccctacccccagaatctctgagtctctgacctaagcaagggatagcaatatggtgagcagttggtagggtttgaagcttaaaatctttaaatgcaaaattgaatttaaatattttaagcttcaaactacaaataccaatttcttttggattcaaggggatcttcatgcagtagacctggttttaagatcagatgatgttgagttattatttcctcattgggagagtagttctaagtaggaggcaggatagcatagctcacggtaccaaaccgatctttagttttttgtttttgtttttcttttgagatgaagtcttactctgtca ccaggctggagtggaatggcgtgcatctcagctcattgcaacctcctcctcccgggttcaagcgattctcccacctca",
          "caaagaagtctggatgggaaatccaaggtattacaattgactccagaaaac aaagggttgaagcacagggtggagggggctgaagccatgtggctcaatggccatgtgaaagaaacaggccaaccAttagtcccagttcactcagtgggAgcctgatacagacggccccacagtgaaagacgggttaacaaaagtcaagagctcaggataagattcttctctggcgtaccaagatttttgttcttggtgccacacctcaccGgatactgac aaagcaaagcaggcaatcc agagaagtgggaaggctcaaacctatgtattaagaagatctggaaagattgggaattttaggcatgaagaagaaagaCtcagaaaggacaagtcccacataaactgccatgtgggaggacatggaacttgtattctatggccactgagattaggttgtaggactaaaggagagaagctatagcaaagggattagggcttaatgagggaggactttcttctgcatgtactgcagcaagAaTgagggagggggatagatgaccacttttcagggagttggaagagattcattcttccaagggatggtcgggcaagacagcctcaaggtcctt cctacccccagaatctctgagtctctgacctaaAcaagggatagcaaCatggtgagcagttggta  gtttgaagcttaaaatctttaaatgcaaaattgaatttaaatattttaagcttcaaactacaaataccaatttctttAggattcaaggggatcttcatgcagtagacctggttttaagatcagatgatgttgagttattatttcctcattgggagagtagttctaagtaggaggcaggatagcatagctca ggtaccaaaccgatctttagttttttgtttttgtttttcttttgagatgaagtcttactctgtcacccaggctggagtggaatggcgtg atctcagctcattgcaacctcctcctcccgggttcaagcgattctcccacctcag",
          20, 20, false, 10);
    final String actStr =         "===========================I=========================================================D================X=======================X===============================================================================================================X=========I===================I=============================================================================X============================================================================D==================================================================================X=X==========================================================================================I=================================X=============X=================II=============================================================================X=================================================================================================================I=================================================================D=======================I======================================================";
    //                       Note: 7 deletes (one of length 2), 3 inserts, 7 SNPs, 1 MNP of length 3 = alignment score 7*2+1 + 3*2 + 7 + 3 = 31.
    assertEquals(actStr, ActionsHelper.toString(actions));
    assertEquals(30, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
    assertEquals(24, actions[ActionsHelper.TEMPLATE_START_INDEX]);
  }

  static final String ACTS1 =          "XX==============D=============I===================================D==============================================X================================================================================D=====================I======================================================I==========================XXX==========================================I================I=======================X====X=============================I==============================================================================III===================================X==================================================================================================II=======================================================X===========================================================================================================================================================================X========================================================================================================X=================================";
  static final String READ1 =          "CAatgtttgagaagta aataggtatcttaccagatgttttgaagtgtatgagaacgatcatctct agccttgacatcagggaagcggccgaacgattatcccccctatattTtctatttcgctctcagaacttggatttcagggtcgcactgccgacatgccgtaccgaattccagatcatcttgcacgact cgcacctgggcacatctgacgccctctaaaaacatctgatgctccaaggtataaaagcaaaggtaggtacttatggttttataccgtttagaaatcgccagaaAATgctatctgctaaggcccaaccacgatatgccgctgagtccgagcctgaccagcatgtcaccattggagtgcagcattaaaaggGaagcCgaatggcagcgaccactaagctggagaggagactaagctgccagccggttcgccctaaacacgtcatggagcaacttgtttgtcgtgaggcggtacttacattctctagggggtattgccgtgataggcaggtatccgatcatcttTagtcgggatagttcggtaccctcagacagactggcgagattgagtgccgctaacaagaaacgtacccaaattgatctcattgtttaaggcagggagccctctctcatgagccttgctgaaatacgctacgttaggtaacaatccttataaagtgtCaatgcacttatgatcaaaatggcattggtcataacgtttgatccacacttgctgaggctggttgttttttcgcgtattgactggctctggaatatcgctgtaggatgcccacaagttgtcattcgtacctagggcacgaactaataatctaatgcccgttgaatatcatgtAaagccatatcggagaggtgcgatggaggtcccacgccgccgcgactcacagcccgtatttggcttagatgactattctgaacgatctcggtttcggctcaggaaGgtccaaacataggtttcatctcacatagtactc";
  static final String TMPL1 = "ctagtgagatgatgtttgagaagtataataggtatctta cagatgttttgaagtgtatgagaacgatcatctctaagccttgacatcagggaagcggccgaacgattatcccccctatattGtctatttcgctctcagaacttggatttcagggtcgcactgccgacatgccgtaccgaattccagatcatcttgcacgactccgcacctgggcacatctgacg cctctaaaaacatctgatgctccaaggtataaaagcaaaggtaggtacttatgg tttataccgtttagaaatcgccagaattggctatctgctaaggcccaaccacgatatgccgctgagtccga cctgaccagcatgtca cattggagtgcagcattaaaaggAaagcGgaatggcagcgaccactaagctggagagg gactaagctgccagccggttcgccctaaacacgtcatggagcaacttgtttgtcgtgaggcggtacttacattctcta   ggtattgccgtgataggcaggtatccgatcatcttGagtcgggatagttcggtaccctcagacagactggcgagattgagtgccgctaacaagaaacgtacccaaattgatctcattgtttaaggcagggagcc  ctctcatgagccttgctgaaatacgctacgttaggtaacaatccttataaagtgtTaatgcacttatgatcaaaatggcattggtcataacgtttgatccacacttgctgaggctggttgttttttcgcgtattgactggctctggaatatcgctgtaggatgcccacaagttgtcattcgtacctagggcacgaactaataatctaatgcccgttgaatatcatgtGaagccatatcggagaggtgcgatggaggtcccacgccgccgcgactcacagcccgtatttggcttagatgactattctgaacgatctcggtttcggctcaggaaAgtccaaacataggtttcatctcacatagtactccgatagtgggcgaaactgg";

  /** Caused AlignmentResult.mismatches() to crash with out of bounds error. */
  public void test1000mismatchesBug() {
    final int[] actions = align(READ1, TMPL1, 10, 26, false, 10);
    assertEquals(ACTS1, ActionsHelper.toString(actions));
    assertEquals(37, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
    assertEquals(9, actions[ActionsHelper.TEMPLATE_START_INDEX]);
  }

  public void testAboveMaxScore() {
    mMaxScore = 36; // override default
    final int[] actions = align(READ1, TMPL1, 10, 0, false, 10);
    assertEquals("", ActionsHelper.toString(actions));
    assertEquals(Integer.MAX_VALUE, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
    assertEquals(9, actions[ActionsHelper.TEMPLATE_START_INDEX]);
  }

  /** Had a lower score than GotohEditDistance. */
  public void test1000TooLow() {
    //                                                   123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789
    final String read =                                "gtgagtacatacagaggaataaagaacactgggaatgcctactatctttgtatctatggcaactatcttgtacatagtaaaatagagatagctgaataaaatgtttttaaataaagcaattCtgataattgacaggcagatatttaagaacCAcagctaaatatttttttctccttagctttgttccatgtcGacattccaaaacaatttcataggcattattttcataattttaatttattttacaattattgttgaatagcaacatttttaattttggacaactaaaacaaatttttcattctattatcctgtttttaactgagttgattttattctttaatatacagacttgccctttcaaaaggcaatacttaaatggaaattttatttccgcatgccagaaccaccaatcacaactgacccggtttacctgttgctcacaagaatctgcaggaaattgaaccctaggccttgcaactcattggactcacaaacactcactgcactggaccccctggggtcaactttatgtcataagtatctgtcaggtatcacagaattaagtgtccagcttataaggatgagtaacattgatttatccctccagaaggaagcttaactttttccttcagtttgctttaatgaagcctcatccatcttgcatgctcaggcagtgagtaatttctgttaagtgaatttttggaagctaaaaattacaaacaaatcttctcattttTttgccttttcttaaaaaaatctaccaaattttccctgacttccaaaattataatgctgtattctagagatgatgtcacttcattgtttcatcaaaagtgagcccctggaatgatgactcatacaggttaatttgcctttaatgtggcaatttgatgcttagaaatggttctgatttattaatgactatctttgtttcccctagaaataaagtctgtgtatagaacaagatccctccatcaatgctcctaaatgttatgga";
    final String tmpl = "agatgggaggaaaatgttgacagatggaaatgtgagtacatacagaggaataaagaacactgggaatgcctactatctttgtatctatggcaactatcttgtacatagtaaaatagagatagctgaataaaatgtttttaaataaagcaattttgataattgacaggca ata ttaagaacatcagctaaatatttttttctccttagctttgttccatgtctacattccaaaacaatttcataggcattattttcataattttaatttattttacaattattgttgaatagcaacatttttaattttggacaactaaaacaaatttttcattctattatcctg ttttaactgagttgattttattctttaatatacagacttg cctttcaaaaggcaatacttaaatggaaattttatttccgcatgccagaaccaccaatcacaactga  cggtttacctgttgctcacaagaatctgcaggaaattgaaccctaggccttgcaactcattggactcacaaacactcactgcactggaccccctggggtcaactttatgtcataagtatctgtcaggt tcacagaattaagtgtccagcttataaggatgagtaacattgatttatccctccagaaggaagcttaactttttccttcagtttgctttaatgaagcctcatccatcttgcatgctcaggcagtgagtaatttctgttaagtgaatttttggaagctaaaaattac aacaaatcttctcattttattgcctt   ttaaaaaaatctaccaaattttccctgacttccaaaattataatgctgtattctagagatgatgtcacttcattgtttcatcaaaagtgagcccctggaatgatgactcatacaggttaatttgcctttaatgtggcaatttgatgcttagaaatggttctgatttattaatgactatctttgtttcccctagaaataaagtctgtgtatagaacaagatccctccatcaatgctcctaaatgttatggacattattgcttttgtctatt";
    //                   0123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 12345678 9 1 23456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 12 3456789 123456789 123456789 123456789 12 3456789 123456789 123456789 123456789 123456789 123456789 123456789   123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 1234567 89 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123 456789 123456789 123456789    123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789
    //                             10        20        30                  50                            80                 100                                               150                                                 200                                                                                                 300                                                                                                   400                                                                                                   500                                                                                                  600                                                                                                 700                                                                                                     800
    final int[] actions = align(read, tmpl, 14, 16, false, 18);
    final String acts =                                "=========================================================================================================================X================I===I========XX=======================================X=========================================================================================================================I========================================I===================================================================II================================================================================================================================I======================================================================================================================================================================I==================X=======III==========================================================================================================================================================================================================================================================";
    assertEquals(acts, ActionsHelper.toString(actions));
    assertEquals(24, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
    assertEquals(31, actions[ActionsHelper.TEMPLATE_START_INDEX]);
    align(read, tmpl, 13, 16, true, 4); // 31 - 13 = 18, which is too far to move

    // now compare with GotohEditDistance
    // (the problem was that the start position was so wrong that Gotoh tried to shift too far)

    final NgsParams params = new NgsParamsBuilder().gapOpenPenalty(1).gapExtendPenalty(1).substitutionPenalty(1).unknownsPenalty(0).create();
    final UnidirectionalEditDistance ed = new GotohEditDistance(params);
    final byte[] readDNA = DnaUtils.encodeString(read.replaceAll(" ", ""));
    final byte[] tmplDNA = DnaUtils.encodeString(tmpl.replaceAll(" ", ""));
    final int[] edSays = ed.calculateEditDistance(readDNA, readDNA.length, tmplDNA, 14, 50, MaxShiftUtils.calculateDefaultMaxShift(readDNA.length), true);
    assertNotNull(edSays);
    //System.out.println("Ed says: " + ActionsHelper.toString(edSays));
    final String edActs = "=========================================================================================================================X================I===I========XX=======================================X=========================================================================================================================I========================================I===================================================================II================================================================================================================================I======================================================================================================================================================================I==================X=======III==========================================================================================================================================================================================================================================================";
    assertEquals(edActs, ActionsHelper.toString(edSays));
    assertEquals(24, edSays[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
    assertEquals(31, edSays[ActionsHelper.TEMPLATE_START_INDEX]);
    final ActionsValidator av = new ActionsValidator(1, 1, 1, 0);
    assertTrue(av.isValid(edSays, readDNA, readDNA.length, tmplDNA, 50));
  }

  public void testRepeatRegion() {
    final String read = "attacaggca tgagccaccg tgccagcctc ttactatttt tttttttttt tttaagatgg ag";
    final String tmpl = "attacaggca tgagccaccg tgccagcctc ttactatttt tttttttttt tttaagatgg aggctcacac tgt";
    align(read, tmpl, 4, 0, false, 4);
    // This is a subtle test: with a shift of -maxShift, the "tttttttttt"
    // in the read uniquely matches the rightmost "tttttttttt" in the template.
    // Which is why findSeed now checks one past maxShift, to ensure that a
    // match is never part of a repeat region.
    align(read, tmpl, MaxShiftUtils.calculateDefaultMaxShift(read.length()) + 1, 0, true, 4);
  }

  // should not handle an indel of length 5.
  public void testIndel5() {
    final String d30 = "aaaaacccccgggggtttttNNaaccggtt";
    align(d30 + d30 + d30 + d30 + "aaaaa" + d30,
          d30 + d30 + d30 + d30 + d30 + d30, 0, 5, true, 4);
  }

  public void testLogStats() {
    final ByteArrayOutputStream bos = new ByteArrayOutputStream();
    Diagnostic.setLogStream(new PrintStream(bos));
    testalign1(); // success
    testalign5(); // cannot initially synchronize
    testalign5(); // cannot initially synchronize
    testalign8(); // cannot resynchronize (10..19) after indel, due to a repeat region
    testalignStartTooHard();
    testalignMiddleTooHard();
    testalignEndTooHard();
    mMaxScore = 50;
    testalign20OffTemplate(); // off template right
    testalign20OffTemplate(); // off template right again
    testalign20OffTemplateLeft();
    testalign31(); // two many indels too close
    testalign24(); // success, but with gotoh call at start and end.
    test1000mismatchesBug(); // success, but with 4 gotoh calls in mid, plus 1 at start.
    testAboveMaxScore();
    testIndel5(); // fails due to indel of length 5.
    mHopStep.logStats();
    Diagnostic.setLogStream();
    final String log = bos.toString();
    //System.out.println(log);
    TestUtils.containsAll(log,
        "HopStepEditDistanceLong[resync_length=7, max_indel=4, indel_spacing=21] statistics:",
        "Total: 15",
        "Success: 4",
        "Cannot synch  : 2",
        "Off template  : 3",
        "Close Indels  : 1",
        "End too hard  : 1 / 2",
        "Mid too hard  : 1 / 7 with 4 merges",
        "Start too hard: 1 / 7",
        "0..: 0",
        "10..: 1",
        "20..: 0",
        "90..: 0");
    assertFalse(log.contains("120..: 1"));
    assertFalse(log.contains("210..: 0"));
    assertFalse(log.contains("220..: 0"));
  }

  public void testExtendSeedRight() {
    final byte[] read = {1, 2, 3, 4, 2, 3, 4, 1, 2, 3};
    final byte[] tmpl = {1, 2, 3, 4, 2, 3, 4, 1, 2, 3, 4, 2, 1};
    //                     indexes: 0  1  2  3  4  5  6  7  8  9 10 11 12
    mHopStep.setupTestData(read, tmpl);
    final SeedPositions seed = new SeedPositions();
    seed.mType = ActionsHelper.SAME;
    seed.mX1 = 0;
    seed.mX2 = 4;
    seed.mY1 = 0;
    seed.mY2 = 4;
    mHopStep.extendSeedRight(seed, 10, 13);
    assertEquals(10, seed.mX2);
    assertEquals(10, seed.mY2);

    seed.mX1 = 1;
    seed.mX2 = 4;
    seed.mY1 = 4;
    seed.mY2 = 7;
    mHopStep.extendSeedRight(seed, 10, 13);
    assertEquals(4, seed.mX2);
    assertEquals(7, seed.mY2);

    seed.mX1 = 0;
    seed.mX2 = 3;
    seed.mY1 = 7;
    seed.mY2 = 10;
    mHopStep.extendSeedRight(seed, 10, 13);
    assertEquals(5, seed.mX2);
    assertEquals(12, seed.mY2);

    // now repeat the previous one, but with a x bound
    seed.mX1 = 0;
    seed.mX2 = 3;
    seed.mY1 = 7;
    seed.mY2 = 10;
    mHopStep.extendSeedRight(seed, 4, 13);
    assertEquals(4, seed.mX2);
    assertEquals(11, seed.mY2);

    // now repeat the previous one, but with a y bound
    seed.mX1 = 0;
    seed.mX2 = 3;
    seed.mY1 = 7;
    seed.mY2 = 10;
    mHopStep.extendSeedRight(seed, 10, 11);
    assertEquals(4, seed.mX2);
    assertEquals(11, seed.mY2);
  }

  public void testRegression() {
    final byte[] t = DnaUtils.encodeString("GTAGACGGATGGATGGATGGACAGGTAGGCAGGTAGACGGATGGTTGGACGGACGGATAGGTAGATGGGTGGGTGGATGGGCGGGTGGATGGGTGGATGGATGGACAGGT");
    final byte[] r = DnaUtils.encodeString("GGATGGATGGATGGTGGGACAGGTAGGCAGGGAGACGGATGGTTGGACGGCGGATAGGTACGATGGGTGGGTGGATGGGCGGGTGGATGGGTGGATGGAT");

    final UnidirectionalEditDistance ed = getEditDistanceInstance(1, 1, 1, 0);
    final int[] actions = ed.calculateEditDistance(r, r.length, t, 2, 15, 7, false);
    assertEquals("X==X==========XX===============X==================D==========I=======================================", ActionsHelper.toString(actions));
  }

//  public void testRegression2() {
//    final byte[] t = DnaUtils.encodeString("GATGGTATCTCATTGTGGTTTTGATTTGCACTTCTCTGATGGCCAGTGATGATGAGCATTTTTTCATGTGTTTTTTGGCTGCATAGATGTCTTCTTTTGAGAAGTGTCTGTTCATATCCT");
//    final byte[] r = DnaUtils.encodeString("             TGTGGTTTTGATTTGCATTTCTCTGATGGCCAGTGATGATGAGCATTTTTTCATGTGTTTTTTGGCTGCATAAATGTCNNNNNNNNNNNNGAAAAAAATG".replaceAll(" ", ""));
//    final byte[] rrev = DnaUtils.encodeString("CATTTTTTTCNNNNNNNNNNNNGACATTTATGCAGCCAAAAAACACATGAAAAAATGCTCATCATCACTGGCCATCAGAGAAATGCAAATCAAAACCACA");
//
//    final HopStepEditDistanceLong ed = new HopStepEditDistanceLong(1, false);
//    final int[] actions = ed.calculateEditDistance(r, r.length, t, 13, 10, 8, false);
//    assertNull(actions);
////
////    System.err.println(Utils.calculateDefaultMaxShift(r.length));
//    final UnidirectionalEditDistance hopl = new HopStepEditDistanceLong(1, false);
//
//    final boolean USELOOPING = false;
//    final EditDistance red = new RcEditDistance(new UnidirectionalPrioritisedEditDistance(hopl, new SeededAligner(false), USELOOPING ? (UnidirectionalEditDistance) new UnidirectionalLoopingEditDistance(new GotohEditDistance(1, false)) : new GotohEditDistance(1, false)));
////    final int[] actions3 = red.calculateEditDistance(r, r.length, t, 13, false, 10, 8, false);
////    System.err.println(ActionsHelper.toString(actions3));
////
////
////    final GotohEditDistance ged = new GotohEditDistance(1, false);
////    final int[] actions2 = ged.calculateEditDistance(r, r.length, t, 13, 10, 8, false);
////    assertEquals(10, actions2[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
////    assertEquals("=================X======================================================X==================XXXXXXX=X", ActionsHelper.toString(actions2));
////
//    System.err.println("rev");
//    final int[] actionsrev = red.calculateEditDistance(rrev, r.length, t, 7, true, 10, 8, false);
//    System.err.println(actionsrev[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
//    System.err.println(ActionsHelper.toString(actionsrev));
//  }

  public void testRcRegression() {

    //                                      0123456789 123456789 1234567 89 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789
    final byte[] t = DnaUtils.encodeString("TAATGCTATCCCTCCCCTAGCCCCCCAC ACTCGATGGGCCCCGGTGTGTGATGTTCCCCTCCCTGTGTCCATGTGTTCTCATTGTTCAACTCCCACTTATGAGTG".replaceAll(" ", ""));
//                                                                                      ------------------- --------------------------------------
//                                            ======================X===IX=X=X=XX======X===================X======================================
    final byte[] r =   DnaUtils.encodeString("ATGCTATCCCTCCCCTAGCCCCACACTCCCCAACAGGCCCCAGTGTGTGATGTTCCCCTCCATGTGTCCATGTGTTCTCATTGTTCAACTCCCACTTATG".replaceAll(" ", ""));
    //                                        0123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789
    //                                                  10        20        30        40        50        60        70        80        90
    // trev = DnaUtils.encodeString("CACTCATAAGTGGGAGTTGAACAATGAGAACACATGGACACAGGGAGGGGAACATCACACACCGGGGCCCATCGAGTGTGGGGGGCTAGGGGAGGGATAGCATTA");
    // rrev                              CATAAGTGGGAGTTGAACAATGAGAACACATGGACACATGGAGGGGAACATCACACACTGGGGCCTGTTGGGGAGTGTGGGGCTAGGGGAGGGATAGCAT

    final int[] actions = mHopStep.calculateEditDistance(r, r.length, t, 2, 10, 8, false);
    assertNull(actions);
    //System.err.println(actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);

    final NgsParams params = new NgsParamsBuilder().gapOpenPenalty(1).gapExtendPenalty(1).substitutionPenalty(1).unknownsPenalty(0).create();
    final GotohEditDistance ged = new GotohEditDistance(params);
    final int[] actions2 = ged.calculateEditDistance(r, r.length, t, 2, 10, 8, false);
    assertEquals(10, actions2[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
    assertEquals(2, actions2[ActionsHelper.TEMPLATE_START_INDEX]);
    assertEquals("======================X===IX=X=X=XX======X===================X======================================", ActionsHelper.toString(actions2));
  }

  public void testFindSeed() {
    final byte[] read = {1, 2, 3, 4, 2, 2, 3, 4, 1, 2, 3};
    final byte[] tmpl = {1, 2, 3, 4, 2, 2, 3, 4, 1, 2, 3, 4, 2, 1};
    //          indexes: 0  1  2  3  4  5  6  7  8  9 10 11 12
    mHopStep.setupTestData(read, tmpl);

    // a shift of 0
    assertEquals("[3 10) [3 10) ", mHopStep.findSeed(6, 6, 9, 9, 4, 1, 0).toString());  //no movement along read

    // some shifts of 0 that fail
    assertNull(mHopStep.findSeed(3, 6, 6, 8, 4, 1, 0));

    // a shift of +2 (with a N at the end)
    assertEquals("[1 5) [9 13) ", mHopStep.findSeed(4, 10, 10, 12, 4, 1, 0).toString());
    assertNull(mHopStep.findSeed(4, 8, 7, 11, 4, 1, 0)); // ambiguous because +3 and -4 both match.

    // a shift of -3
    assertNull(mHopStep.findSeed(3, 7, 9, 10, 4, 1, 0)); // ambiguous because -3 and +4 both match.
    assertNull(mHopStep.findSeed(3, 6, 9, 9, 4, 1, 0));

    // a shift of +2 with an N
    assertEquals("[1 5) [9 13) ", mHopStep.findSeed(4, 10, 7, 13, 4, 1, 0).toString());

    // a shift of 0
    assertEquals("[0 5) [8 13) ", mHopStep.findSeed(5, 12, 9, 13, 4, 2, 2).toString()); // move down at most once

    // now check that some entropy is needed to find a match.
    final byte[] boring = {1, 2, 3, 4, 4, 4, 4, 4};
    //                       indexes: 0  1  2  3  4  5  6  7
    mHopStep.setupTestData(boring, boring);
    assertNull(mHopStep.findSeed(6, 6, 6, 9, 4, 1, 0)); // no movement along read
    assertEquals("[2 7) [2 7) ", mHopStep.findSeed(6, 6, 6, 9, 4, 2, 1).toString());
  }

  public void testDifferent() {
    assertTrue(mHopStep.different((byte) 1, (byte) 4));
    assertFalse(mHopStep.different((byte) 1, (byte) 1));
    assertTrue(mHopStep.different((byte) 0, (byte) 1));
    assertTrue(mHopStep.different((byte) 1, (byte) 0));
    assertTrue(mHopStep.different((byte) 0, (byte) 0));
  }

  public void testPenaltiesMismatchTrue() {
    final UnidirectionalEditDistance ed = getEditDistanceInstance(19, 1, 9, 0);

    final byte[] r = DnaUtils.encodeString("AtttttgttttCtgccccccccc");
    final byte[] t = DnaUtils.encodeString("gtttttgttttgtgccccccccc");

    final int[] actions = ed.calculateEditDistance(r, r.length, t, 0, 20, 5, false);
    assertNotNull(actions);
    assertEquals(18, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
    assertEquals(0, actions[ActionsHelper.TEMPLATE_START_INDEX]);
    assertEquals("X==========X===========", ActionsHelper.toString(actions));
  }
}
