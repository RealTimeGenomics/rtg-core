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

import com.rtg.mode.DnaUtils;
import com.rtg.ngs.NgsParams;
import com.rtg.ngs.NgsParamsBuilder;
import com.rtg.util.MaxShiftUtils;
import com.rtg.util.diagnostic.Diagnostic;

import junit.framework.TestCase;

/**
 */
public class LowerBoundEstimatorTest extends TestCase {

  private LowerBoundEstimator align(final String readstring, final String templatestring, final int zeroBasedPos, final int score, final int maxScore, int unknownsPenalty) {
    Diagnostic.setLogStream();
    final byte[] read = DnaUtils.encodeString(readstring);
    final byte[] template = DnaUtils.encodeString(templatestring);
    final LowerBoundEstimator ed = new LowerBoundEstimator(1, unknownsPenalty);
    final int lb = ed.calcLB(read, read.length, template, zeroBasedPos, maxScore, MaxShiftUtils.calculateDefaultMaxShift(read.length));
    final int gs = gotohscore(readstring, templatestring, zeroBasedPos, unknownsPenalty);
    //System.err.println("gotoh:" + gs + " lb:" + lb);

    assertEquals(score, lb);

    assertTrue(gs + " >= " + lb, gs >= lb);
    return ed;
  }

  public void test4eq() {
    final String s1 = "ctag";
    final String s2 = "ctag";
    align(s1, s2, 0, 0, 0, 1);
  }

  public void test4eqn1() {
    final String s1 = "ctnn";
    final String s2 = "ctag";
    align(s1, s2, 0, 0, 0, 1);
  }

  public void test4eqn2() {
    final String s1 = "ctag";
    final String s2 = "ctnn";
    align(s1, s2, 0, 0, 0, 1);
  }

  public void test4eqm1() {
    final String s1 = "nnag";
    final String s2 = "ctag";
    align(s1, s2, 0, 0, 0, 1);
  }

  //bug1353
  public void test4eqm2() {
    final String s1 = "ctag";
    final String s2 = "ntag";
    align(s1, s2, 0, 0, 0, 1);
  }

  public void test8eq() {
    final String s1 = "ctagatga";
    final String s2 = "ctagatga";
    align(s1, s2, 0, 0, 0, 1);
  }

  public void test4eq1x3eq() {
    final String s1 = "ctagatga";
    final String s2 = "ctagacga";
    final LowerBoundEstimator ed = align(s1, s2, 0, 1, 1, 1);
    assertEquals(0, ed.getCount(0));
    assertEquals(1, ed.getCount(1));
    assertEquals(1, ed.getCount(2));
    assertEquals(1, ed.getCount(3));

  }

  public void test4eq2x2q() {
    final String s1 = "ctagatga";
    final String s2 = "ctagtcga";
    align(s1, s2, 0, 1, 1, 1);
  }

  public void test4eq2x2q2() {
    final String s1 = "ctagataa";
    final String s2 = "ctagtcga";
    align(s1, s2, 0, 1, 1, 1);
  }

  public void test4eq2x2q2MI() {
    final String s1 = "ctagataa";
    final String s2 = "ctagtcga";
    align(s1, s2, 0, 1, 0, 1);
  }

  public void test4eq2x2q3() {
    final String s1 = "ctagataa";
    final String s2 = "ctagtcgg";
    align(s1, s2, 0, slow(s1, s2, 1), 1, 1);
  }

  public void test4eq2x2q4() {
    final String s1 = "ctagataat";
    final String s2 = "ctagtcggg";
    final LowerBoundEstimator ed = align(s1, s2, 0, slow(s1, s2, 2), 2, 1);
    assertEquals(2, ed.getCount(0));
    assertEquals(1, ed.getCount(1));
    assertEquals(1, ed.getCount(2));
    assertEquals(1, ed.getCount(3));
  }

  public void test4eq2x2q4x() {
    final String s1 = "ctagatgatttt";
    final String s2 = "ctagtcgacccc";
    align(s1, s2, 0, slow(s1, s2, 2), 2, 1);
  }

  public void test4eq2x2q4xmix() {
    final String s1 = "ctagatgattttaataaata";
    final String s2 = "ctagtcgaccccaataaata";
    align(s1, s2, 0, slow(s1, s2, 3), 3, 1);
  }

  public void test4eq2x2q4xmixMV() {
    final String s1 = "ctagatgattttaataaata";
    final String s2 = "ctagtcgaccccaataaata";
    align(s1, s2, 0, slow(s1, s2, 2), 2, 1); // max gap of 2 (window size)
  }

  public void test4eq2x2q4xmixScore4() {
    final String s1 = "ctagatgatNtttaataaata";
    final String s2 = "ctagatgattttaataaata";
    align(s1, s2, 0, 0, 4, 1);
  }

  public void test4eq2x2q4xmixScore4b() {
    final String s1 = "ctagatgatttttaataaata";
    final String s2 = "ctagatgatttttaataaata";
    align(s1, s2, 0, slow(s1, s2, 0), 0, 1);
  }

  public void test3eq1x() {
    final String s1 = "ctag";
    final String s2 = "ctac";
    align(s1, s2, 0, slow(s1, s2, 1), 1, 1);
  }

  public void test3eq1x1eq() {
    final String s1 = "ctaga";
    final String s2 = "ctaca";
    align(s1, s2, 0, slow(s1, s2, 1), 1, 1);
  }

  public void test3eq1x2eq() {
    final String s1 = "ctagat";
    final String s2 = "ctacat";
    align(s1, s2, 0, slow(s1, s2, 1), 1, 1);
  }

  public void test3eq1x3eq() {
    final String s1 = "ctagatt";
    final String s2 = "ctacatt";
    align(s1, s2, 0, slow(s1, s2, 1), 1, 1);
  }

  public void test3eq1x4eq() {
    final String s1 = "ctagatta";
    final String s2 = "ctacatta";
    align(s1, s2, 0, slow(s1, s2, 1), 1, 1);
  }

  public void test3eq1x4eq1x() {
    final String s1 = "ctagattan";
    final String s2 = "ctacattat";
    align(s1, s2, 0, 0, 3, 1);
  }

  public void test3eq1x4eq1xnswapped() {
    final String s1 = "ctagattat";
    final String s2 = "ctacattan";
    align(s1, s2, 0, 0, 1, 1);
  }

  public void test3eq1x4eq1xlotstemplateNs() {
    final String s1 = "ctagattattactgatgg";
    final String s2 = "ctacattannnnnnnnnn";
    align(s1, s2, 0, 0, 1, 1);
  }

  public void test3eq1x4eq1xlotsreadNs() {
    final String s1 = "ctagattatnnnnnnnnnnn";
    final String s2 = "ctacattattactgatgg";
    align(s1, s2, 0, 0, 1, 1);
  }

  public void testMaxScore1() {
    final String s1 = "ctaga";
    final String s2 = "ctaca";
    align(s1, s2, 0, slow(s1, s2, 1), 1, 1);
  }

  public void testMaxScore1Offset5() {
    final String s1 = "ctaga";
    final String s2 = "nnnnnctaca";
    final LowerBoundEstimator ed = align(s1, s2, 5, 1, 1, 1);
    assertEquals(1, ed.getCount(0));
    assertEquals(0, ed.getCount(1));
    assertEquals(0, ed.getCount(2));
    assertEquals(1, ed.getCount(3));
  }

  public void testMaxScore1Offset50() {
    final String s1 = "ctaga";
    final String s2 = "nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnctaca";
    align(s1, s2, 50, 1, 0, 1);
  }

  public void testMaxScore1Offset50perfect() {
    final String s1 = "ctaga";
    final String s2 = "nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnctaga";
    align(s1, s2, 50, 0, 99, 1);
  }

  public void testMaxScore1Offset50perfectalln() {
    final String s1 = "ctaga";
    final String s2 = "nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn";
    align(s1, s2, 50, 0, 0, 1);
  }

  public void testMaxScore1Offset50perfectalln2() {
    final String s1 = "nnnnn";
    final String s2 = "nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn";
    align(s1, s2, 50, 0, 0, 1);
  }

  public void testMaxScore1Offset50perfectalln2neg50() {
    final String s1 = "nnnnn";
    final String s2 = "nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn";
    align(s1, s2, -50, 0, 0, 1);
  }

  public void testMaxScore1Offset50perfectalln2neg50b() {
    final String s1 = "nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn";
    final String s2 = "nnnnnnnnnnnnnnnnnnnnnnn";
    align(s1, s2, -50, 0, 0, 1);
  }


  public void testMaxScore1b() {
    final String s1 = "ttata";
    final String s2 = "ctatta";
    align(s1, s2, 0, slow(s1, s2, 1), 1, 1);
  }

  public void testReal1() {
    final String s1 = "TTCGTATCTCGATCGGCCGACTACGGGCTCTAGGGCTACCGAAATCCCGCCATCGCAACTACCTGCTTCAATTTTGCAGGACCTACCGACCTGCTAGCGG";
    final String s2 = "tgactgatcgatcgatcgatcgatcgatcgatgctagctagctagctagcTTCGTATGTCGATCGGCCGACCACGGGCTCTTGGGCTACCGAAATCCCGCCATCGCAACTACCTGCTTCAATTTTGCAGGAGCTACCGACCTGCTAGCGG";
    //                 .12345678901234567890123456789012345678901234567890123456789012345678901234567890
    align(s1, s2, 50, 4, 9999, 1);
  }

  private int slow(final String s1, final String s2, final int maxScore) {
    final int[] count = new int[4];
    for (int i = 0; i < s1.length() - 3; ++i) {
      final String ts = s1.substring(i, i + 4);
      //      System.err.println(ts);
      final int low = Math.max(0, i - maxScore);
      final int high = Math.min(i + maxScore + 4, s2.length());
      final String tem = s2.substring(low, high);
      if (!tem.contains(ts)) {
        //       System.err.println("can't find: " + ts + " in " + tem);
        count[i % 4]++;
      }
    }
    int lb = 0;
    for (final int element : count) {
      //      System.err.println(i+" "+count[i]);
      if (element > lb) {
        lb = element;
      }
    }
    return lb;
  }

  public void testReal2() {
    final String s1 = "GGCTCTCTACGATAGGTATTTGTGCGGGTGTAGCTCCACGTATGCTGCTCACACCCAACGCGCCGGGCCCCCCCTCGTTCTAGCATTATACCGGGGTTGCCCAAAGTCTCAGGGTCGACTACTGGCTGTAGGCATCGTACGA CTTTTCCGCTTTGTTCTCTGCGAGAGGCCTGATCGAAAACAGACTCTTTTAATAATTTTACTCTGCCAATAATGTGCAGAGGGGTGTTGCAGTCGAATTGATGAGAGCAAATCTTACTGTAGTATTAGTAGGCCTGTCGTT GGTCATAGCAACCCGGATCCTACTCTAGGCCTTCGACGGAAGGATCACTT TCCCTTTGCGAGCAATGACCAATGGGTAACTGTCGTAGGGACCGTGTCCCAGAGTCTGGTTGTCGGCTACGATCCTCATGTATCACATGAGCAATTGCTACAACGGCCCGGCCACGGTTGTAGGGCGGTTCGTGTAACCCGCCATATAGGCACATTTCGAATAGACAGAGTTCCATTGTCTAAAGTTCGGAGCGCAAAAATACGGTTGCATTCGTGATTGCTGCCCGTGCAAGGCTCGACTACTCCTCCAACGACCGTGAGCCAAGGCGACTCACGCTGCACAACTCGCTGATATACTACACAGATACCGTAAGCAGACACCCCATGGTTCCCACCAATCGGATCACAATATGCACGGGGGCCTCTAACACGCTGTAGTCGGTAGCATCCATTAGTTTAAGTAGTCGAGCCGTTAGCCTAATCTCGTTGGAGAGAGAATCAGCATTTTCAGATCGTTGGTATTCTCGTTTTGGGACTAACAATTTCAGTTATTCCGTTCACCTTCCAAGCGTTCCAGCTAGACCGTCCGTCTGTGGGCATAATCTCTCACTTTATGAAAGCGGTGTAACTCAACTGAGGGGCGGATACCTGGTATCAGGACACCCCTGCCGGGGTTGCAGTCTACTGGTCGCGGCCTTGTGACGGTTTGCGGGAACCAGAATTTAA";
    final String s2 = "GGCTCTCTACGATAGGTATTTGTGCGGGTGTAGCTCCACGTATGCTGCTCACACCCAAACGCGCCGGGCCCCCCTCGTTCTAGCATTATTCTGGGGTTGCC AAAGTCTCAGGGTCGACTACTGGCTGTAGGCATCGTACGATCTTTTCCGCTTTGTTCTCTGCGAGAGGCCTGATCGAAAACAGACTCTTTTAATAATTTTACTCTGCCAATAATGTGCAGAGGGATGTTGCAGTCGAATTGATGAGAGCAAATCTTACTGTAGTATTAGTAGGCCTGTCGTTTGGTGATAGCAACCCGGATCCTACTCTAGGCCTTCGACGGAAGGATACACTTTCCCTTTGCGAGCAATGACCAATGGGTAACTGTCGTAGGGACCGTGTCCCAGAGTCTGGTTGTCGGCTACGATCCTCATGTATCACATGAGCAATTGCTACAACGGCCCGGCCACGGATGTAGGGCGGTTCGTGTAACCCGCCATTAGGCACATTTCGAATAGACAGAGTTCCATTGTCTAAAGTTCGGACGCAAAAATACGGTTGCATTCGTGATTGCTGCCCGTGCAAGGCTCGACTACTCCTCCAACGACCGTGAGCCAAGGCGACTCACGCTGCAGAACTCGCTGATATACTACACAGATACCGTAAGCAGACACCCCATGGTTCACCAATCGGATCACAATATGCACGGGGGCCTCTAACACGCTGTAGTCGGTAGCATCCATTAGTTTAAGTAGTCGAGCCGTTAGCCTAATCTCGTTGGAGAGAGAATCAGCATTTTCAGTTCGTTGGTATTCTCCTTTTGGGACTAACAATCTCAGTTATTCCGTTCTCCTTCCAAGCGTTCCAGCTAGACCGTCCGTCTGTGGGCATAATCTCTCACTTTATGAAAGCGGTGTAACTCAACTGAGGGGCGGATACCTGCTATCAGGACACCCCTGCCGGGGTTGCAGTCTACGGGTCGCGGCCTTGTGACGGTTTGCGGGAACCATAATTTAA";
    //                 .12345678901234567890123456789012345678901234567890123456789012345678901234567890
    //TODO this needs checking by hand - this is a regression test
    align(s1.replaceAll(" ", ""), s2.replaceAll(" ", ""), 0, 18, 100, 1);
  }

  public int gotohscore(final String s1, final String s2, final int zerostart, int unknownsPenalty) {
    final NgsParams params = new NgsParamsBuilder().gapOpenPenalty(1).gapExtendPenalty(1).substitutionPenalty(1).unknownsPenalty(unknownsPenalty).create();
    final UnidirectionalEditDistance ed = new GotohEditDistance(params.gapOpenPenalty(), params.gapExtendPenalty(), params.substitutionPenalty(), params.unknownsPenalty(), false);
    final byte[] read = DnaUtils.encodeString(s1);
    final byte[] template = DnaUtils.encodeString(s2);
    final int[] lb = ed.calculateEditDistance(read, read.length, template, zerostart, 9999, MaxShiftUtils.calculateDefaultMaxShift(read.length), true);
    return lb[ActionsHelper.ALIGNMENT_SCORE_INDEX];
  }



  public void testReal3() {
    final String s1 = "CTCGCGATCCCACCAATCATTCAGACCCCGCAGCAAAGCCCACAGACTTACTCGTACCTGCTGGTAAGAAATCGTAGCACTACGTCCTCGGGAGTCAAAG";
    final String s2 = "CTCGCGATCCCACCAATCATTCAGACCCCGCAAGCACAGACTTACTCGTACCTGCTGGTAAGAATCGTAGCACTGCGTCCTCGGGAGTCAAAGATGTCAG";
    //                 .12345678901234567890123456789012345678901234567890123456789012345678901234567890
    //    slow(s1, s2, 150);
    align(s1, s2, 0, slow(s1, s2, 10), 10, 1);
  }

  public void testReal4vnegative() {
    final String s1 = "CTCGCGATCCCACCAATCATTCAGACCCCGCAGCAAAGCCCACAGACTTACTCGTACCTGCTGGTAAGAAATCGTAGCACTACGTCCTCGGGAGTCAAAG";
    final String s2 = "CTCGCGATCCCACCAATCATTCAGACCCCGCAAGCACAGACTTACTCGTACCTGCTGGTAAGAATCGTAGCACTGCGTCCTCGGGAGTCAAAGATGTCAG";
    //                 .12345678901234567890123456789012345678901234567890123456789012345678901234567890
    //    slow(s1, s2, 150);
    align(s1, s2, -930, 4, 10, 1);
  }

  public void testSimple1() {
    final String s1 = "catgcgtagcg";
    final String s2 = "catgtgcgtagcg";
    //                 0123012301230
    //    slow(s1, s2, 150);
    final LowerBoundEstimator ed = align(s1, s2, 0, 1, 2, 1);
    assertEquals(0, ed.getCount(3));
    assertEquals(1, ed.getCount(0));
    assertEquals(0, ed.getCount(1));
    assertEquals(0, ed.getCount(2));
  }

  public void testSimple2() {
    final String s1 =   "catgcgtagcg";
    final String s2 = "catgtgcgtagcgtgatgctgatcgatc";
    //                 0123012301230
    //    slow(s1, s2, 150);
    final LowerBoundEstimator ed = align(s1, s2, 0, 2, 1, 1);
    assertEquals(2, ed.getCount(0));
    assertEquals(1, ed.getCount(1));
    assertEquals(1, ed.getCount(2));
    assertEquals(1, ed.getCount(3));
  }

  public void testSimple3() {
    final String s1 = "catgcgtagcg";
    final String s2 = "catgtgcgtagcgtgatgctgatcgatc";
    //                 0123012301230
    //    slow(s1, s2, 150);
    final LowerBoundEstimator ed = align(s1, s2, 0, 1, 0, 1);
    assertEquals(1, ed.getCount(0));
    assertEquals(0, ed.getCount(1));
    assertEquals(0, ed.getCount(2));
    assertEquals(0, ed.getCount(3));
  }

  public void testSimple4() {
    final String s1 = "tcat";
    final String s2 = "catg";
    //                 0123012301230
    //    slow(s1, s2, 150);
    final LowerBoundEstimator ed = align(s1, s2, 0, 1, 0, 1);
    assertEquals(0, ed.getCount(0));
    assertEquals(0, ed.getCount(1));
    assertEquals(0, ed.getCount(2));
    assertEquals(1, ed.getCount(3));
  }

  public void testSimple5() {
    final String s1 = "aaaatttt";
    final String s2 = "ttttaaaa";
    //                 0123012301230
    //    slow(s1, s2, 150);
    final LowerBoundEstimator ed = align(s1, s2, 0, 1, 4, 1);
    assertEquals(1, ed.getCount(0));
    assertEquals(1, ed.getCount(1));
    assertEquals(1, ed.getCount(2));
    assertEquals(0, ed.getCount(3));
  }

  public void testSimple6() {
    final String s1 = "aaaatttt";
    final String s2 = "ttttaaaa";
    //                 0123012301230
    //    slow(s1, s2, 150);
    final LowerBoundEstimator ed = align(s1, s2, 0, 1, 3, 1);
    assertEquals(1, ed.getCount(0));
    assertEquals(1, ed.getCount(1));
    assertEquals(1, ed.getCount(2));
    assertEquals(1, ed.getCount(3));
  }

  public void testSimple7() {
    final String s1 = "aaaatttg";
    final String s2 = "ttttaaaa";
    //                 0123012301230
    //    slow(s1, s2, 150);
    final LowerBoundEstimator ed = align(s1, s2, 0, 1, 4, 1);
    assertEquals(1, ed.getCount(0));
    assertEquals(1, ed.getCount(1));
    assertEquals(1, ed.getCount(2));
    assertEquals(1, ed.getCount(3));
  }

  public void testSimple8() {
    final String s1 = "aaaatntn";
    final String s2 = "ttttaaaa";
    //                 0123012301230
    //    slow(s1, s2, 150);
    final LowerBoundEstimator ed = align(s1, s2, 0, 0, 8, 1);
    assertEquals(1, ed.getCount(0));
    assertEquals(0, ed.getCount(1));
    assertEquals(0, ed.getCount(2));
    assertEquals(0, ed.getCount(3));
  }

  public void testSimple9() {
    final String s1 = "aaaatntn";
    final String s2 = "ttttaaaa";
    //                 0123012301230
    //    slow(s1, s2, 150);
    final LowerBoundEstimator ed = align(s1, s2, 0, 0, 1, 1);
    assertEquals(1, ed.getCount(0));
    assertEquals(0, ed.getCount(1));
    assertEquals(0, ed.getCount(2));
    assertEquals(1, ed.getCount(3));
  }

  public void testSimple10() {
    final String s1 = "aaaatntn";
    final String s2 = "ttntaaaa";
    //                 0123012301230
    //    slow(s1, s2, 150);
    final LowerBoundEstimator ed = align(s1, s2, 0, 0, 1, 1);
    //assertEquals(1, ed.getCount(0)); //TODO
    assertEquals(0, ed.getCount(1));
    assertEquals(0, ed.getCount(2));
    //assertEquals(1, ed.getCount(3)); //TODO
  }

  public void testSimpletata() {
    final String s1 = "tatataat";
    final String s2 = "tataaatata";
    //                 0123012301230
    //    slow(s1, s2, 150);
    final LowerBoundEstimator ed = align(s1, s2, 0, 1, 3, 1);
    assertEquals(1, ed.getCount(0));
    assertEquals(0, ed.getCount(1));
    assertEquals(0, ed.getCount(2));
    assertEquals(1, ed.getCount(3));
  }

  public void testSimpletata4() {
    final String s1 = "tatataat";
    final String s2 = "tataaatata";
    //                 0123012301230
    //    slow(s1, s2, 150);
    final LowerBoundEstimator ed = align(s1, s2, 0, 1, 4, 1);
    assertEquals(0, ed.getCount(0));
    assertEquals(0, ed.getCount(1));
    assertEquals(0, ed.getCount(2));
    assertEquals(1, ed.getCount(3));
  }

  public void testSimpletata4b() {
    final String s1 = "tatataat";
    final String s2 = "tataatata";
    //                 0123012301230
    //    slow(s1, s2, 150);
    final LowerBoundEstimator ed = align(s1, s2, 0, 0, 3, 1);
    assertEquals(0, ed.getCount(0));
    assertEquals(0, ed.getCount(1));
    assertEquals(0, ed.getCount(2));
    assertEquals(0, ed.getCount(3));
  }

  public void testSimpleaann1() {
    final String s1 = "aaaaaa";
    final String s2 = "nnnnnnnnnnnnaaaannnnnnnnnnnnnnnn";
    //                 0123012301230
    //    slow(s1, s2, 150);
    final LowerBoundEstimator ed = align(s1, s2, 0, 0, 99, 1);
    assertEquals(0, ed.getCount(0));
    assertEquals(0, ed.getCount(1));
    assertEquals(0, ed.getCount(2));
    assertEquals(0, ed.getCount(3));
  }

  public void testSimpleaann2() {
    final String s1 = "aaataaa";
    final String s2 = "nnnnnnnnnnnnaaaannnnnnnnnnnnnnnn";
    //                 0123012301230
    //    slow(s1, s2, 150);
    final LowerBoundEstimator ed = align(s1, s2, 0, 0, 0, 1);
    assertEquals(0, ed.getCount(0));
    assertEquals(0, ed.getCount(1));
    assertEquals(0, ed.getCount(2));
    assertEquals(0, ed.getCount(3));
  }

  public void testSimpleaann3() {
    final String s1 = "aaaataaa"; //aaaataaa
    final String s2 = "nnngnnnnnnnnaaaanaaannnnnnnnnnnn";
    //                 0123012301230
    //    slow(s1, s2, 150);
    final LowerBoundEstimator ed = align(s1, s2, 0, 0, 99, 1);
    assertEquals(0, ed.getCount(0));
    assertEquals(0, ed.getCount(1));
    assertEquals(0, ed.getCount(2));
    assertEquals(0, ed.getCount(3));
  }

  public void testSimpleaann4() {
    final String s1 = "aaaataaa";
    //                  aaaataaa
    final String s2 = "nnaaaanaaannn";
    //                 0123012301230
    //    slow(s1, s2, 150);
    final LowerBoundEstimator ed = align(s1, s2, 2, 0, 99, 1);
    assertEquals(0, ed.getCount(0));
    assertEquals(0, ed.getCount(1));
    assertEquals(0, ed.getCount(2));
    assertEquals(0, ed.getCount(3));
  }


  public void testSimpleaann4off1() {
    final String s1 = "aaaataaa";
    //                  aaaataaa
    final String s2 = "naaaanaaannn";
    //                 0123012301230
    //    slow(s1, s2, 150);
    final LowerBoundEstimator ed = align(s1, s2, 0, 0, 99, 1);
    assertEquals(0, ed.getCount(0));
    assertEquals(0, ed.getCount(1));
    assertEquals(0, ed.getCount(2));
    assertEquals(0, ed.getCount(3));
  }

  public void testSimpleaann42() {
    final String s1 = "aaaataattannna";
    //                  aaaataaa
    final String s2 = "ttttaaatgaaannn";
    //                 0123012301230
    //    slow(s1, s2, 150);
    final LowerBoundEstimator ed = align(s1, s2, 0, 2, 1, 1);
    assertEquals(1, ed.getCount(0));
    assertEquals(1, ed.getCount(1));
    assertEquals(1, ed.getCount(2));
    assertEquals(2, ed.getCount(3));
  }


  public void testEdges() {
    final String s1 = "gtgtgtgtgtgt";
    final String s2 = "acacacacacaccacacaccacacacaccaca";
    //                 01234567890123456789012345678901
    //                 0         1         2         3
    align(s1, s2, 0, 0, 100, 0);
    align(s1, s2, 1, 0, 100, 0);
    align(s1, s2, 3, 0, 100, 0);
    align(s1, s2, 6, 0, 100, 0); //maxshift within front end
    align(s1, s2, 7, 3, 100, 0);
    align(s1, s2, 12, 3, 100, 0);
    align(s1, s2, 13, 3, 100, 0);
    align(s1, s2, 14, 0, 100, 0); //maxshift within rear end
    align(s1, s2, 17, 0, 100, 0);
    align(s1, s2, 20, 0, 100, 0);
  }

}
