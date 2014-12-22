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
import com.rtg.util.PortableRandom;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.MemoryPrintStream;

/**
 */
public class SeededAlignerTest extends AbstractUnidirectionalEditDistanceTest {

  @Override
  protected UnidirectionalEditDistance getEditDistanceInstance(int gapOpenPenalty, int gapExtendPenalty, int subsPenalty, int unknownsPenalty) {
    final NgsParams params = new NgsParamsBuilder().gapOpenPenalty(gapOpenPenalty).gapExtendPenalty(gapExtendPenalty).substitutionPenalty(subsPenalty).unknownsPenalty(unknownsPenalty).create();
    return new SeededAligner(params, false);
  }

  public int[] align(String readstring, String templatestring, int zeroBasedPos, int score, boolean nullisok, int maxScore, Integer expectedStart, int gapOpenPenalty, boolean treatNsAsMismatches) {
    final byte[] read = DnaUtils.encodeString(readstring.replaceAll(" ", ""));
    final byte[] template = DnaUtils.encodeString(templatestring.replaceAll(" ", ""));
    final UnidirectionalEditDistance ed = getEditDistanceInstance(gapOpenPenalty, 1, 1, treatNsAsMismatches ? 1 : 0);
    final int[] actions = ed.calculateEditDistance(read, read.length, template, zeroBasedPos, maxScore, MaxShiftUtils.calculateDefaultMaxShift(read.length), true);

    if (!nullisok) {
      assertNotNull(readstring + " " + templatestring, actions);
    }
//        System.err.println(actions[ActionsHelper.TEMPLATE_START_INDEX]);

    if (actions != null) {
//       System.err.println(actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
//       System.err.println("1: " + readstring + "\n2: " + templatestring);
      assertEquals(readstring + " " + templatestring, score, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
      if (expectedStart != null) {
        assertEquals(expectedStart.intValue(), actions[ActionsHelper.TEMPLATE_START_INDEX]);
      }
    }
    return actions;
  }

  /**
   * To test alignment speed, set repeatCount above to 100,000 then run this
   * test.
   */
  public void testNoSeeds() {
    align("atgtacgtacgtacgtacgtacgtacgtacgtacgt", "acgtacgtacgNacgtacgtacgtacgtacgtacgt", 0, 0, true, 100, null, 1, false);
  }

  public void testOneSeedMiddle() { //01234567890123456789012345678901234567890
    align("atcggctagtcatgagcta",     "nnnnnnnnnnnnnnnnnnnnatcggctagtcatgagctannnnnnnnnnnnnnnnnnnnnn", 24, 0, false, 100, null, 1, false);
  }

  public void testOneBigContiguousSeed() {
    //012345678901234567890123456789012345678901234567890
    align("atcggctagtcattgtagctagtcttcgcgcttagatgatagaaacgtta", "atcggctagtcattgtagctagtcttcgcgcttagatgatagaaacgtta", 0, 0, false, 100, null, 1, false);
  }

  public void testOneSeed() {
    //012345678901234567890123456789012345678901234567890
    align("atcggctagtcattgtagctagtcttcgcgcttagatgatagaaacgtta", "nnatcggctagtcattgtagctagtcttcgcgcttagatgatagaaacgttann", 0, 0, false, 100, null, 1, false);
  }

  public void testTwoSeeds() {
    //012345678901234567890123456789012345678901234567890
    align("atcggctagtcattgtagctagtcttcgcgcttagatgatagaaacgtta", "nnnnatcggctagtcattgtagctagtctnnntcgcgcttagatgatagaaacgttannnnnn", 0, 4, false, 100, null, 1, false);
  }

  public void testTwoSeedsFewNs() {
         //0         1         2         3         4         5
         //012345678901234567890123456789012345678901234567890       atcggctagtcattgtagc  tagtct   tcgcgcttagatgatagaaa  cgtta
    align("atcggctagtcattgtagctagtcttcgcgcttagatgatagaaacgtta", "nnnnatcggctagtcattgtagcnntagtctnnntcgcgcttagatgatagaaanncgttannnnnn", 4, 10, false,
        100, null, 1, false);
  }

  public void testTwoSeedsFewNs2() {
    //012345678901234567890123456789012345678901234567890            atcggctagtcattgtagc  tagtct   tcgcgcttagatgatagaaa  cgtta
    align("atcggctagtcattgtagctagtcttcgcgcttagatgatagaaacgtta", "nnnnatcggctagtcattgtagcnntagtctnnntcgcgcttagatgatagaaanncgttannnnnn", 4, 10, false,
        10, null, 1, false);
  }

  public void testTwoSeedsFewNs2Termearly() {
    //012345678901234567890123456789012345678901234567890
    align("atcggctagtcattgtagctagtcttcgcgcttagatgatagaaacgtta", "nnnnatcggctagtcattgtagcnntagtctnnntcgcgcttagatgatagaaanncgttannnnnn", 4,
        Integer.MAX_VALUE, false, 9, null, 1, false);
  }

  public void testTwoSeedsFewNs2MaxScore2() {
    //012345678901234567890123456789012345678901234567890
    align("atcggctagtcattgtagctagtcttcgcgcttagatgatagaaacgtta",
      "nnnnatcggctagtcattgtagcnntagtctnnntcgcgcttagatgatagaaanncgttannnnnn", 4,
        Integer.MAX_VALUE, true, 2, null, 1, false);
  }


  // a case that intermittently fails on .NET...
  // TODO Come back to this test and actually make it a valid test.
  // The seeded aligner returns NULL if there is 10nt repeated regions in the two strings.
  // Hence the null check fails.
//  public void testVariousSidesSameKnownFail() {
//    for (int s = 0; s < 100; s++) {
//      final PortableRandom r = new PortableRandom(s);
//      final char[] dna = {'a', 't', 'c', 'g'};
//      final long start = System.currentTimeMillis();
//      int count = 0;
//      while (System.currentTimeMillis() - start < 100) {
//        for (int i = SeededAligner.NTHASHSIZE; i < 100; i += 3) {
//          final StringBuilder sb = new StringBuilder();
//          for (int j = 0; j < i; j++) {
//            sb.append(dna[r.nextInt(dna.length)]);
//          }
//          align(sb.toString(), sb.toString(), 0, 0, i < SeededAligner.NTHASHSIZE ? true : false, Integer.MAX_VALUE);
//          count++;
//        }
//      }
//    }
//  }


  // 100ms of random strings with no Ns against the same string with randomly inserted Ns
  public void testVariousSidesSameWithNsAgainstAnotherStringwithNs() {
    final PortableRandom r = new PortableRandom();
    final char[] dna = {'a', 'a', 'a', 'a', 't', 't', 't', 't', 't', 'c', 'c', 'c', 'c', 'c', 'g', 'g', 'g', 'g'};
    final long start = System.currentTimeMillis();
    while (System.currentTimeMillis() - start < 100) {
      for (int i = SeededAligner.NTHASHSIZE; i < 100; i += 3) {
        final StringBuilder sb = new StringBuilder();
        for (int j = 0; j < i; j++) {
          sb.append(dna[r.nextInt(dna.length)]);
        }
        final StringBuilder sb2 = new StringBuilder(sb.toString());
        for (int j = 0; j < sb2.length(); j++) {
          if (r.nextInt(5) < 1) {
            sb2.setCharAt(j, 'n');
          }
        }
        //       System.err.println(sb+" "+sb2);
        align(sb.toString(), sb.toString(), 0, 0, true, Integer.MAX_VALUE, null, 1, false);
      }
    }
  }

  // 100ms of random strings with no Ns against the same string with 1 random sub ONLY
  public void testInsertingOneRandomSub() {
    final PortableRandom r = new PortableRandom();
    final char[] dna = {'a', 't', 'c', 'g'};
    final long start = System.currentTimeMillis();
    while (System.currentTimeMillis() - start < 100) {
      for (int i = SeededAligner.NTHASHSIZE + 1; i < 60; i += 3) {
        final StringBuilder sb = new StringBuilder();
        for (int j = 0; j < i; j++) {
          sb.append(dna[r.nextInt(dna.length)]);
        }
        final StringBuilder sb2 = new StringBuilder(sb.toString());
        int pos = 0;
        while (pos < SeededAligner.NTHASHSIZE) {
            pos = r.nextInt(sb2.length());
        }
        final char ch = sb2.charAt(pos);
        char ch2 = ch;
        if (i > SeededAligner.NTHASHSIZE) {
          if (ch == 'a') {
            ch2 = 't';
          } else if (ch == 't') {
            ch2 = 'a';
          } else if (ch == 'c') {
            ch2 = 'g';
          } else if (ch == 'g') {
            ch2 = 'c';
          }
        }
        sb2.setCharAt(pos, ch2);
        //       System.err.println(sb+" "+sb2);
        align(sb.toString(), sb2.toString(), 0, 1, true, Integer.MAX_VALUE, null, 1, false);
      }
    }
  }


  public void testExact15() {
    //    .12345678901234567890123456789012345678901234567890
    align("atcggctagtcattg", "atcggctagtcattgt", 0, 0, false, 100, null, 1, false);
  }

  public void testExact14() {
    //    .12345678901234567890123456789012345678901234567890
    align("atcggctagtcatt", "atcggctagtcattgt", 0, 0, true, 100, null, 1, false);
  }

  public void testExact14temp0() {
    //    .12345678901234567890123456789012345678901234567890
    align("atcggctagtcatt", "", 0, 0, true, 100, null, 1, false);
  }

  public void testAllNs() {
    //    .12345678901234567890123456789012345678901234567890
    align("nnnnnnnnnnnnnnnnnnnn", "nnnnnnnnnnnnnnnnnnnn", 0, 0, true, 100, null, 1, false);
  }

  public void testbigread2nsreplacedeitherendinthetemplate() {
    //    .12345678901234567890123456789012345678901234567890
    align("atcggctagtcattatgagtacga", "nncggctagtcattatgagtacnn", 0, 0, false, 100, null, 1, false);
  }

  public void testbigread2missingoffeitherendinthetemplate() {
    //    .12345678901234567890123456789012345678901234567890
    align("atcggctagtcattatgagtacga", "cggctagtcattatgagtac", 0, 4, false, 100, null, 1, true);
  }

  public void testtemplate2eitherendlonger() {
    //    .12345678901234567890123456789012345678901234567890
    align("atcggctagtcattatgagtacga", "taatcggctagtcattatgagtacgaat", 0, 0, false, 100, null, 1, false);
  }

  public void testtwoseeds1subgap() {
    //    .12345678901234567890123456789012345678901234567890
    align("atcggctagtcattatgagtcacgattgagatgctatgttta", "atcggctagtcattatgagttacgattgagatgctatgttta", 0, 1, false, 100, null, 1, false);
  }

  public void testtwoseeds1delinread() {
    //    .12345678901234567890123456789012345678901234567890
    align("atcggctagtcattatgagtacgattgagatgctatgttta", "atcggctagtcattatgagttacgattgagatgctatgttta", 0, 2, false, 100, null, 1, false);
  }

  public void testtwoseeds1insinread() {
    //    .12345678901234567890123456789012345678901234567890
    align("atcggctagtcattatgagtttacgattgagatgctatgttta", "atcggctagtcattatgagttacgattgagatgctatgttta", 0, 2, false, 100, null, 1, false);
  }

  public void testoneseedoffsetby5() {
    //    .12345678901234567890123456789012345678901234567890
    align("atcggctagtcattatgagt", "ttattatcggctagtcattatgagt", 0, 0, false, 100, null, 1, false);
  }

  public void testoneseedoffsetbynegative5() {
    //    .12345678901234567890123456789012345678901234567890
    align("atcggctagtcattatgagt", "ctagtcattatgagtatgatga", 0, 5, false, 100, -5, 1, true);
  }     //                    atcggctagtcattatgagt

  public void testoneseedoffset5start5() {
    //    .12345678901234567890123456789012345678901234567890
    align("atcggctagtcattatgagt", "ttattatcggctagtcattatgagt", 5, 0, false, 100, null, 1, false);
  }

  public void testoneseedoffset5start6() {
    //    .12345678901234567890123456789012345678901234567890
    align("atcggctagtcatta", "tttttatcggctagtcatta", 6, 0, true, 100, null, 1, false);
  }

  public void testoneseedoffset5start5b() {
    //    .12345678901234567890123456789012345678901234567890
    align("atcggctagtcatta", "tttttatcggctagtcatta", 5, 0, false, 100, null, 1, false);
  }

  public void testfourseedswithNs() {
    //    .12345678901234567890123456789012345678901234567890123456789012345678901234567890
    align("atcggctagtcattannnatgcatgctattagannngctggatgctgacgannntgatgtgatatatggc",
          "atcggctagtcattaatgcatgctattagagctggatgctgacgatgatgtgatatatggc", 0, 3 * (3 + 1), false, 100, null, 1, false);
  }

  public void testfourseedswithNsInTemplate() {
    //    .12345678901234567890123456789012345678901234567890123456789012345678901234567890
    align("atcggctagtcattaatgcatgctattagagctggatgctgacgatgatgtgatatatggc",
        "atcggctagtcattannnnnatgcatgctattagannnnngctggatgctgacgannnnntgatgtgatatatggc", 0, 3 * (5 + 1), false, 100, null, 1, false);
  }

  public void testfourseedswithNsInTemplateMax5() {
    //    .12345678901234567890123456789012345678901234567890123456789012345678901234567890
    align("atcggctagtcattaatgcatgctattagagctggatgctgacgatgatgtgatatatggc",
        "atcggctagtcattannnnnatgcatgctattagannnnngctggatgctgacgannnnntgatgtgatatatggc", 0, Integer.MAX_VALUE, false, 5, null, 1, false);
  }

  public void testfourseedswithNsInTemplateMax4() { // tight maxscore, so tight we don't even find seeds
    //    .12345678901234567890123456789012345678901234567890123456789012345678901234567890
    align("atcggctagtcattaatgcatgctattagagctggatgctgacgatgatgtgatatatggc",
        "atcggctagtcattannnnnatgcatgctattagannnnngctggatgctgacgannnnntgatgtgatatatggc", 0, Integer.MAX_VALUE, true, 4, null, 1, false);
  }

  public void testfourseedswithNsInTemplateMax18() {
    //    .12345678901234567890123456789012345678901234567890123456789012345678901234567890
    align("atcggctagtcattaatgcatgctattagagctggatgctgacgatgatgtgatatatggc",
        "atcggctagtcattannnnnatgcatgctattagannnnngctggatgctgacgannnnntgatgtgatatatggc", 0, 18, false, 18, null, 1, false);
  }

  public void testfourseedswithNsInTemplateMax17() {
    //    .12345678901234567890123456789012345678901234567890123456789012345678901234567890
    align("atcggctagtcattaatgcatgctattagagctggatgctgacgatgatgtgatatatggc",
        "atcggctagtcattannnnnatgcatgctattagannnnngctggatgctgacgannnnntgatgtgatatatggc", 0, Integer.MAX_VALUE, false, 17, null, 1, false);
  }

  //  2010-02-25 12:14:28            TGGGGATCCAGGACCCTTTAGCAACTACTAGCCTTTTTTTTTTTTTTTTTTTTTTTTTCTTTTTTGCCAGATTCTTCTCCTTAGCAACTGTAATGGCATA
  //  2010-02-25 12:14:28              AAATTCATCTTGCATAGTCCACAAAGTTCTGTCATTTTTTTTTTTTTTTTTTTTGAGATGGAGTCTCACTCTGTCACCGAGGCTGCCGGGCAGTGGCACC

  public void testlargereal1() {
    //    .12345678901234567890123456789012345678901234567890123456789012345678901234567890
       align("TGGGGATCCAGGACCCTTTAGCAACTACTAGCCTTTTTTTTTTTTTTTTTTTTTTTTTCTTTTTTGCCAGATTCTTCTCCTTAGCAACTGTAATGGCATA",
        "AAATTCATCTTGCATAGTCCACAAAGTTCTGTCATTTTTTTTTTTTTTTTTTTTGAGATGGAGTCTCACTCTGTCACCGAGGCTGCCGGGCAGTGGCACC", 0, Integer.MAX_VALUE, true, 100, null, 1, false);
  }

  public void testlargereal2() {
    final String s1 = "AGACTCCATCTCAAAATAAATAAATAAATAAATAAATAAATAAATAAATAAATAAATATTTTATGATCAAAGATACTGTCCTTGGGAAGCTAAATTGGAC";
    final String s2 = "CAAGACTCTGTCTCAATAAATAAATAAATAAATAAACAAACAAACAAACAGTTTAACATAAATATATGTATTACTTGGCTGGTAGGATTACAGGTGATTG";
    //    .12345678901234567890123456789012345678901234567890123456789012345678901234567890
//    slowCheck(s1, s2);
    align(s1, s2, 0, Integer.MAX_VALUE, true, 30, null, 1, false);
  }

  public void testlargereal3() {
    final String s1 = "TGGGGATCCAGGACCCTTTAGCAACTACTAGCCTTTTTTTTTTTTTTTTTTTTTTTTTCTTTTTTGCCAGATTCTTCTCCTTAGCAACTGTAATGGCATA";
    final String s2 = "CTCATTAGCATGTTAAATATCCACCCAGGTTTTTTTTTTTTTTTTTTTTTTTTTACTATTATCATCAGCAAAGGGTCAGCTTGAGGGCAGGTTAAATCAA";
    //    .12345678901234567890123456789012345678901234567890123456789012345678901234567890
//    slowCheck(s1, s2);
    align(s1, s2, 0, Integer.MAX_VALUE, true, 100, null, 1, false);
  }

// This test is not high enough complexity :(
//  public void testlargereal4() {
//    final String s1 = "aactacccgcccttccgatgatagtgttggtatgccctgcgcagaacgatgaaaagttgggtacgtgtcccagtagcagctttccggtttccattatgagactttccattggggagtcgcagccctaaacttgtttaccttggctacccgaaagtcggatgattggatgagagcaaggctggtaggcctcaacgtagtcttttacaaatgatacggaagggcattcatgtaggtcatcacgctgtcaatggaccgtcacatcttttat";
//    final String s2 = "aactacccgcccttccgatgatagtgttggtatgccctgcgcagaacgatgaaaagttgggtacgtgtcccagtagcagctttccggtttccattatgagactttccataggggagtcgcagccctaaacttgtttaccttggctacccgaaagtcggatgattggatgagagcaaggctggtaggcctcaacgtagtcttttacaaaagatacggaagggcattcatgtaggtcatcacgctgtcaatggaccgtcacatcttttat";
//    //    .12345678901234567890123456789012345678901234567890123456789012345678901234567890
////    slowCheck(s1, s2);
//    align(s1, s2, 0, 2, false, Integer.MAX_VALUE, null);
//  }

  public void testlargereal6() {
    final String s1 = "ATGTCCCTACAAAGGACATGAACTCATCATTTTTTANGGCTGCATAGTATTCCATGGTGCATATGTGACAC";
    final String s2 = "ATGTCCCTACAAAGGACATGAACTCATCATTTTTTATGGCTGCATAGTATTCCATGGTGTATATGTGCCAC";
    //                 .12345678901234567890123456789012345678901234567890123456789012345678901234567890
    //                                      123456789012345
    // slowCheck(s1, s2);
    // checkSubs(s1, s2);
    align(s1, s2, 0, 2, false, Integer.MAX_VALUE, null, 1, false);
  }

  public void testlargereal7() {
    final String s1 =  "CACCTATGAGTGAGAACATGCGGTGTTTGGTTTTTTGTCCTTGTGATAGTTTGCTAAGAATGATGGTTTCCAGCTTCATCCATGTCCCTACAAAGGACATGAACTCATCATTTTTTANGGCTGCATAGTATTCCATG";
    final String s2 = "CCACCTATGAGTGAGAACATGCGGTGTTTGGTTTTTGTCCTTGCAATAGTTTGCTGAGAATGATGGTTTCCAGTTTAATCCATGTCCCTACAAAGGACATGAACTCATCATTTTTTATGGCTGCATAGTATTCCATG";
    //    .12345678901234567890123456789012345678901234567890123456789012345678901234567890
    // slowCheck(s1, s2);
    // checkSubs(s1, s2);
    align(s1, s2, 0, 7, false, Integer.MAX_VALUE, null, 1, false);
  }


  public void testlargereal8() {
    final String s1 = "ATGTCCCTACAAAGGACATGAACTCATCATTTTTTANGGCTGCATAGTATTCCATGGTGCATATGTGACAC";
    final String s2 = "ATGTCCCTACAAAGGACATGAACTCATCATTTTTTATGGCTGCATAGTATTCCATGGTGTATATATGCCAC";
    //                 .12345678901234567890123456789012345678901234567890123456789012345678901234567890
    // slowCheck(s1, s2);
    // checkSubs(s1, s2);
    align(s1, s2, 0, 3, false, Integer.MAX_VALUE, null, 1, false);
  }


  public void testlargereal9() {
    final String s1 = "AATGTGGCACATATACACCATGGAATACTATGCAGCCNTAAAAATGATGAGTTCTATGTCCTTTGTAGGGACATGGATGAAGCTGGAAAC";
    final String s2 = "AATGTGGCACATATACATCATGGAATACTATGCAGCCATAAAAATGATGAGTTCATGTCCTTTGTAGGGACATGGATGAAATTGGAAATC";
    //    .12345678901234567890123456789012345678901234567890123456789012345678901234567890
    // slowCheck(s1, s2);
    // checkSubs(s1, s2);
    align(s1, s2, 0, 6, false, Integer.MAX_VALUE, null, 1, false);
  }

  public void testFixedBothMaxScore() {
    final String s1 = "AATGTGGCACATATACATCATGGAATACTATTATATTTATTTTTATTATTAATCATGTCCTTTGTAGGGACATGGATGAAATTGGAAATC";
    final String s2 = "AATGTGGCACATATACATCATGGAATACTATGCAGCCATAAAAATGATGAGTTCATGTCCTTTGTAGGGACATGGATGAAATTGGAAATC";
    //    .12345678901234567890123456789012345678901234567890123456789012345678901234567890
    // slowCheck(s1, s2);
    // checkSubs(s1, s2);
    align(s1, s2, 0, Integer.MAX_VALUE, true, 6, null, 1, false);
  }
  public void testObviousDelete() {
    final String s1 = "AATGTGGCACATATACATCATGGAATACATTTGTAGGGACATGGATGAAATTGGAAATC";
    final String s2 = "AATGTGGCACATATACATCATGGAATACTTTGTAGGGACATGGATGAAATTGGAAATC";
    //    .12345678901234567890123456789012345678901234567890123456789012345678901234567890
    // slowCheck(s1, s2);
    // checkSubs(s1, s2);
    align(s1, s2, 0, 2, true, 6, null, 1, false);
  }

  public void testlargereal10() {


//                                            ACTTTAAGTTTTAGGG
//                                               TTAAGTTTTAGGGTACATGTGCACA
//                                                                         TGTGCAGGTTAGTTACATATGT
//                                                                                    GTTACATATGTATAC
//                                                                                       ACATATGTATACATGTGCCATGCTGGTG
    final String s1 = "TTTTCTTTTTTTATTATTT   TACTTTAAGTTTTAGGGTACATGTGCACAGTGTGCAGGTTAGTTACATATGTATACATGTGCCATGCTGGTGCGCTGCACC".replaceAll(" ", "");
    final String s2 = "TTTTCTTTTTTTATTATTATTATACTTTAAGTTTTAGGGTACATGTGCACAATGTGCAGGTTAGTTACATATGTATACATGTGCCATGCTGGTGTGCTGCACC".replaceAll(" ", "");
//                                            ACTTTAAGTTTTAGGG
//                                               TTAAGTTTTAGGGTACATGTGCACA
//                                                                         TGTGCAGGTTAGTTACATATGT
//                                                                                    GTTACATATGTATAC
//                                                                                       ACATATGTATACATGTGCCATGCTGGTG
//                                           TACTTTAAGTTTTAGGGTACATGTGCACA TGTGCAGGTTAGTTACATATGTATACATGTGCCATGCTGGTG
    //                 .123456789112345678921234567893123456789412345678951234567896123456789712345678981234567899123456789012
//    slowCheck(s1, s2);
    // checkSubs(s1, s2);
    align(s1, s2, 0, 6, false, Integer.MAX_VALUE, null, 1, false);
    //gotohcheck(s1, s2, 0, 6);
    //boundedcheck(s1, s2, 0, 6);
  }


  public void testlargereal12() {
    final String s1 = removeSpaces("GGGAATGTGGTTCATCTAG GGAGCGT GGCCCATCTAGGGAG CGTGGTTCATCTAGGGAGCGT GGCCCATCTAGGGAG CGTGGTTCATCTAGGGAGCGTGG");
    final String s2 = removeSpaces("GGGAGTGTGGTTCATCTA GGGAGCGTGGCCCATCTAG GGAATGTGGTTCATCTAGGGAGCGT  GGCCCATCTAGGGAG GGTGGCTCATCTAGGGCGCGTGG");
    //                 .12345678901234567890123456789012345678901234567890123456789012345678901234567890
//    slowCheck(s1, s2);
    // checkSubs(s1, s2);
    align(s1, s2, 0, 6, true, 20, null, 1, false);
//    gotohcheck(s1, s2, 0, 6);
//    boundedcheck(s1, s2, 0, 6);
  }

  public void testNumCoveredTooLow() {
    final String s1 = removeSpaces("ATCTAG GGAGCGT GGCCCATCAATGCGT GGCCCATCTAGGGAG C");
    final String s2 = removeSpaces("GGGAGTGTGGTTCATCTA GGGAGCGTGGCCCATCTAG GGAATGTGGTTCATCTAGGGAGCGT  GGCCCATCTAGGGAG GGTGGCTCATCTAGGGCGCGTGG");
    align(s1, s2, 6, -1, true, 100, null, 1, false);
  }


  public void testMaxShiftWeirdness() {
    final String read = "CCCTTGTAAGTTGGATTCCTAGGTATTTTATTCTCTTTGAAACAATTGTGAATAGAAGTTCACTCATGATTTGGCTCTCTGTTTGTCTGTTGTTGGTGNA";
    final String tmpl = "CCCTTGTAAGTTGGATTCCTAGGTATTTTATTTTCTTTGAAGCAATTGTGAATGGGAGTTCACTCATGATTTGGCTGTTTGTCTGTTATTGGTGTATAAG";

    align(read, tmpl, 0, 10, false, 100, 0, 1, false);
  }


  String removeSpaces(String s) {
    final StringBuilder sb = new StringBuilder();
    for (int i = 0; i < s.length(); i++) {
      if (s.charAt(i) != ' ') {
        sb.append(s.charAt(i));
      }
    }
    return sb.toString();
  }

  public void testlargereal11() {
    // Both Gotoh and CachedMatrix edit distances give a global score of 16, not 17.
    // The seeded aligner+gotoh gives 16, but seed aligner+cachedmatrix gives 17.
    //                                v                  v                                                      v   v                v v       v
    final String s1 = removeSpaces("CATTCACAATTGCTTCAAA     TAA AATACCTAGGAATCCAACTTACAAGGGATGTGAAGGACCTCTTCAAG GAGAA CTACAAACCACTGCTCAGGGAAATAG");
    final String s2 = removeSpaces("ACAATTGCTTCAAAGAGAATAA AATACCTAGGAATCCAACTTACAAGGGATGTGAAGGACCTCTTCAAG CAGAG CTACAAACCACTGCTGAAGGAAATAA");
    //                 .12345678901234567890123456789012345678901234567890123456789012345678901234567890
//    slowCheck(s1, s2);
    // checkSubs(s1, s2);
    align(s1, s2, 0, 16, false, Integer.MAX_VALUE, null, 1, true);
    //gotohcheck(s1, s2, 0, 16);
    //boundedcheck(s1, s2, 0, 17);
  }

  public void testlargereal13() {
    // Gotoh     gives 22 (Cigar: 7X10=1X5=1X28=7I  1=2X6=1I3=1X34=),
    // Cach.Mat. gives 23 (Cigar: 7X10=1X5=1X27=7I1X1=2X6=1I3=1X34=).
    //                                v8             v9    v10                            v13       v       v23
    final String s1 = removeSpaces("GTCTCAGATACAAAATAAATGTG CAAAAATCACAAGCATTCTTATACACCA ATTAACAGA        AGC C AAATCATGAGTGAACTCCCATTCACAATTGCTTC");
    final String s2 = removeSpaces("      GATACAAAATCAATGTA CAAAAATCACAAGCATTCTTATACACCA AC AACAGACAAGA G AGC A AAATCATGAGTGAACTCCCATTCACAATTGCTTC");
//    slowCheck(s1, s2);
    // checkSubs(s1, s2);
    align(s1, s2, 0, 19, false, Integer.MAX_VALUE, null, 1, true);
    //gotohcheck(s1, s2, 0, 22);
    //boundedcheck(s1, s2, 0, 17);
  }

  public void testSeedPositions() {
    final SeedPositions pos = new SeedPositions();
    assertEquals("[0 0) [0 0) ", pos.toString());
  }

  public void testIntermittentFail() {
    final String s1 = "cgttaatcccgcgtggggatctgcatggataagctgaccggtctgtagtctccccccccccccaaacacagagtttctctaatcgccgcacact";
    final String s2 = "cgttaatcccgcgtggggatctgcatggataagctgaccggtctgtagtctccccccccccccaaacacagagtttctctaatcgccgcacact";
    align(s1, s2, -1, 1, true, 1, null, 1, false);
  }

  public void testMaxShift() {
    try (MemoryPrintStream mps = new MemoryPrintStream()) {
      Diagnostic.setLogStream(mps.printStream());                                                                                                                                                                                                         //     atcggctagtcattgtagctagtcttcgcgcttagatgatagaaacgtta                                                                                                                                                                          attcatgtaggtcatcacgctgtcaatggacc
      final String tmpl = "AGTAAGGGGG ctccctgccccctgggggaacctgaccctaattctgactccagtagcagctttccagtagcagctttccagtagcagctttccagtagcagctttccagtagcagctttccagtagcagctttccagtagcagctttccagtagcagctttgcaggggg aagcaggcgctcgaggggctcaca gatgcacacacagtagcagctttccagtagcagcatcggctagtcattgtagctagtcttcgcgcttagatgatagaaacgttaagcagctttccagtagcagctttccagtagcagctttcttcacagattat acctattagtggtcaagacaaacacgtggctctcagatctaaaccttctttagggacactaaagaaactg aggcttgggacagacagcagtagcagctttccagtagcagctttccagattcatgtaggtcatcacgctgtcaatggaccagctttccagtagcagctttccagtagcagctttccagtagcagcttttgaggccaggacaccaagccagccctgggccctcactgctcagccaggacccttttctcctgaactaaaccatctgggttt agttga actaaacctggggctgcccgtgcaaaggcatctgccacattcagttctcagcaagttggcaacactgacaagcatccctgcaagagaggtgaagactcaagttctatttgtttctggcatg aaaagcgatggttt   ccttctcatcactcataaaaaatagcagaaatgtggcaattgcgggatcgatcccaccgcacaccatacgtgtgcatatgctaccatgcacataggcctcaacccaggcatgggaaataccagtaggagcggcaaagagaccttccaatctctctagaaactccctacaactgaaccagtgctacttcagatcccgccgcacaccatacgcgtgcatatgctaccatgcacataggacgtgaggtctggcgtcaactcaggcatgggaaataccagtaggagtggcaaagagaccttccaatctctctagaaactccctacaactggaccagcgctacttcaggtgacagccgcgtacattcctcaaatgatgacaacaggcagccagaataccgc".replaceAll(" ", "");
      final String read = "agtaaggggggctccctgccccctgggggaacctgaccctaattctgactccagtagcagctttccagtagcagctttccagtagcagctttccagtagcagctttccagtagcagctttccagtagcagctttccagtagcagctttccagtagcagctttgcaggggggaagcaggcgctcgaggggctcacaagatgcacacacagtagcagctttccagtagcagcatcggctagtcattgtagctagtcttcgcgcttagatgatagaaacgttaagcagctttccagtagcagctttccagtagcagctttcttcacagaatcgaacctattagtggtcaagacaaacacgtggctctcagatctaaaccttctttagggacactaaagaaactggaggcttgggacagacagcagtagcagctttccagtagcagctttccagattcatgtaggtcatcacgctgtcaatggaccagctttccagtagcagctttccagtagcagctttccagtagcagcttttgaggccaggacaccaagccagccctgggccctcactgctcagccaggacccttttctcctgaactaaaccatctgggttttagttgacactaaacgtggggctgcccgtgcaaaggcatctgccacattcagttctcagcaagttggcaacactgacaagcatccctgcaagagaggtgaagactcaagttctatttgtttctggcatggaaaagcgatggttttttccttctcatcactcataaaaaatagcagaaatgtggcaattgcgggatcgatcccaccgcacaccatacgtgtgcatatgctaccatgcacataggcctcaacccaggcatgggaaataccagtaggagcggcaaagagaccttccaatctctctagaaacttcctacaactgaaccagtgctacttcagatcccgccgcacaccatacgcgtgcatatgctaccatgcacataggacgtga".replaceAll(" ", "");
      //          -                                                                                                                                                               -                        -                                                                                                                                  -----                                                                      -                                                                                                                                                                                                                 -      -                                                                                                                         -              ---                                                                                                                                                                  -
      final String readShort = "ggctccctgccccctgggggaacctgaccctaattctga";

      align(readShort, tmpl, 0, 99, true, 10, 8, 1, false);
      align(readShort, tmpl, 1, 0, false, 10, 8, 1, false);
      align(readShort, tmpl, 7, 0, false, 10, 8, 1, false);
      align(readShort, tmpl, 10, 0, false, 10, 8, 1, false);
      align(readShort, tmpl, 15, 0, false, 10, 8, 1, false);
      align(readShort, tmpl, 16, -1, true, 10, 8, 1, false);

      align(tmpl, tmpl, -1, -1, true, 32, 0, 1, false);

      final UnidirectionalEditDistance ed = getEditDistanceInstance(1, 1, 1, 1);
      int[] actions = ed.calculateEditDistance(DnaUtils.encodeString(tmpl), tmpl.length(), DnaUtils.encodeString(tmpl), 0, 32, MaxShiftUtils.calculateDefaultMaxShift(tmpl.length()), true);
      assertNotNull(actions);
      assertEquals(0, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
      assertEquals(0, actions[ActionsHelper.TEMPLATE_START_INDEX]);
      //UnidirectionalEditDistance ed = align(tmpl, tmpl, 0, 0, false, 32, 0);
      ed.logStats();
      align(tmpl, tmpl, 9, 0, true, 32, 0, 1, false);
      align(tmpl, tmpl, 17, 0, false, 32, 0, 1, false);
      align(tmpl, tmpl, 18, 0, false, 32, 0, 1, false);
      align(tmpl, tmpl, 19, -1, true, 32, 0, 1, false);
      align(tmpl, tmpl, 20, -1, true, 32, 0, 1, false);

      align(read, tmpl, -1, -1, true, 32, 0, 1, false);
      align(read, tmpl, 0, 25, false, 32, 0, 1, false);
      align(read, tmpl, 1, 25, false, 32, 0, 1, false);
      align(read, tmpl, 5, 25, false, 32, 0, 1, false);
      align(read, tmpl, 10, 25, false, 32, 0, 1, false);
      actions = ed.calculateEditDistance(DnaUtils.encodeString(read), read.length(), DnaUtils.encodeString(tmpl), 11, 32, MaxShiftUtils.calculateDefaultMaxShift(read.length()), true);
      assertNull(actions);
      //= align(read, tmpl, 11, 39, true, 32, 0);
      ed.logStats();

      TestUtils.containsAll(mps.toString(), "calls: 2",
        "37      0      0      0      0      0      0      0      0      0",
        "0      3      0      0      0      0      0      0      0      0",
        "0      0      1      0      0      0      0      0      0      0");
    }
  }

  public void testMaxShiftFixedBoth() {
          //             GCCCAGTTAACTAGTCATTGTAGCTAGTCTTACTGGGTTTCACCA            TCCTGACCTCAGGTGAGTCACCTGCCTTGGCCTCC
    final String read = "GCCCAGTTAACTAGTCATTGTAGCTAGTCTTACTGGGTTTCACCACGTTGGCCAGGCTCCTGACCTCAGGTGAGTCACCTGCCTTGGCCTCCCAAAGTGC";
    final String tmpl = "GCCCAGTTAACTAGTCATTGTAGCTAGTCTTACTGGGTTTCACCATGTTGGCCAGGCTGGGGAATTCCTGACCTCAGGTGAGTCACCTGCCTTGGCCTCCCAAAGTGCTGGGGTTACAAGCATGAGTCACCGCACCTGGCCACTTCCCACAG";
    align(read, tmpl, 0, 10, false, 50, 0, 1, false);
  }

  public void testMaxShiftFixedStart() {
    //                   CGCCATCTTGTTTGTGTGCTAGGCTGGGGGGGAGAGAGGGCGAGAGAGAGCGGGCGAGAGTGGGCAAGCAGGACGCCGGGCTG
    final String read = "CGCCATCTTGTTTGTGTGCTAGGCTGGGGGGGAGAGAGGGCGAGAGAGAGCGGGCGAGAGTGGGCAAGCAGGACGCCGGGCTGCTGCGGGAGCCAGAGAG";
    final String tmpl = "CGCCATCTTGTTTGTGTGCTAGGCTGGGGGGGAGAGAGGGCGAGAGAGAGCGGGCGAGAGTGGGCAAGCAGGACGCCGGGCTGAGTGCTAACTGCGGGAGCCAGAGAGTGCGGAGGGGAGTCGGGTCGGAGAGAGGCGGCAGGGGCCGAGACAGTGGCAGGG";
    align(read, tmpl, 0, 9, false, 50, 0, 1, false);

  }

  public void testMNPAlignment() {
    final String read1 =     "ATATCCAAGAAAATATTACCAAATCCAGTGTCATAAAGGTATTCCCCTATGTTTTCTTCTAAGAGTTTTATGTTTTATGTTTAGGTCTTTGACCCATTTCAAGATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGAGACGGAGTCTCGCTCTGTCGCCCAGGTCGGACTGCGGACTGCAGTGGCGCAATCCCGGCTCAC";
    final String template1 = "ATATCCAAGAAAATATTACCAAATCCAGTGTCATAAAGGTATTCCCCTATGTTTTCTTCTAAGAGTTTTATGTTTTATGTTTAGGTCTTTGACCCATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGAGACGGAGTCTCGCTCTGTCGCCCAGGTCGGACTGCGGACTGCAGTGGCGCAATCCCGGCTCAC";

    align(read1, template1, 0, 5, false, 50, 0, 1, false);



//    596     83      gi|89161203|ref|NC_000022.9|NC_000022   139     255     95=8I18=8D79=   =       2284146 -206    ATATCCAAGAAAATATTACCAAATCCAGTGTCATAAAGGTATTCCCCTATGTTTTCTTCTAAGAGTTTTATGTTTTATGTTTAGGTCTTTGACCCATTTCAAGATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGAGACGGAGTCTCGCTCTGTCGCCCAGGTCGGACTGCGGACTGCAGTGGCGCAATCCCGGCTCAC        ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++        AS:i:18 NM:i:16 MQ:i:255        XA:i:37 IH:i:1
//    224561  99      gi|89161203|ref|NC_000022.9|NC_000022   147     255     91=5X99=1X4=    =       2284178 218     GAAAATATTACCAAATCCAGTGTCATAAAGGTATTCCCCTATGTTTTCTTCTAAGAGTTTTATGTTTTATGTTTAGGTCTTTGACCCATTTCAAGATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGAGACGGAGTCTCGCTCTGTCGCCCAGGTCGGACTGCGGACTGCAGTGGCGCAATCCCGGCTCACTGCTAGCT        ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++        AS:i:6  NM:i:6  MQ:i:255        XA:i:24 IH:i:1

    final String read =             "ATATCCAAGAAAATATTACCAAATCCAGTGTCATAAAGGTATTCCCCTATGTTTTCTTCTAAGAGTTTTATGTTTTATGTTTAGGTCTTTGACCCATTTCAAGATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGAGACGGAGTCTCGCTCTGTCGCCCAGGTCGGACTGCGGACTGCAGTGGCGCAATCCCGGCTCAC";
    final String template = "TTGGTGTCATATCCAAGAAAATATTACCAAATCCAGTGTCATAAAGGTATTCCCCTATGTTTTCTTCTAAGAGTTTTATGTTTTATGTTTAGGTCTTTGACCCATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGAGACGGAGTCTCGCTCTGTCGCCCAGGTCGGACTGCGGACTGCAGTGGCGCAATCCCGGCTCACTGCAAGCTCCGCTTCCC";

    align(read, template, 8, 5, false, 50, 8, 1, false);

    align(DnaUtils.reverseComplement(read), DnaUtils.reverseComplement(template), 17, 5, true, 50, 17, 1, false);
  }

  public void testRepetitive() {
    final String read = "ACACATACATACATACACACACACACACACACACACACACACACACACACAAGCCCACATATATATATAT  TTTTTTTTTTTAAGGCAGAGTCTCACTCTG".replaceAll(" ", "");
    final String tmpl = "ACACATACATACATACACACACACACACACACACACACACACACACACA      CACATATATATATATATTTTTTTTTTTTAAGGCAGAGTCTCACTCTGTCTTAAGGCTGGAGTGCAGTGTCGTGATCTCG".replaceAll(" ", "");

    align(read, tmpl, 0, -1, true, 50, 0, 1, false);
  }

  public void testGapOpenPenalty() {
    //one insert, one substitution, 2 gap open penalty
    align("atcggctagtcattatgag tacgattgagatactatgttta",
          "atcggctagtcattatgagttacgattgagatgctatgttta", 0, 4, false, 100, null, 2, false);

    //one delete, one substitution, 2 gap open penalty
    align("atcggctagtcattatgagtttacgattgagatactatgttta",
          "atcggctagtcattatgagtt acgattgagatgctatgttta", 0, 4, false, 100, null, 2, false);

    //one insert, one substitution, 0 gap open penalty
    align("atcggctagtcattatgag tacgattgagatactatgttta",
        "atcggctagtcattatgagttacgattgagatgctatgttta", 0, 2, false, 100, null, 0, false);
  }

  public void testNegSeedLen() {
    final String read = "                              AGGGAGTCTTTTGTCCAATAAACCCAATTTGAAGGCCGGGCCGCGGTGGCTCACACCTGTAATCCCAGCACTATTGGGAGGCCGAGTGGACCCCCCGGATCCATGAGTCAGGAGATCGAGAACCATCCTGGCTAACACGGTGAAACCCTGTTCTCTACCTAAAAAATACAAAAATTAGGCCGGGCGAGGTGGCCGGCCGGCCGTGTAGTCCCAGCTACTTGGGAGGCTGAGGCAGGAGAATGGCAATGAACCCTGGGGGACGGAGCCTGCAGTGAGCTGAGATGCCGCCACTGTACTCCAGCCTGGTGACAGTGAGAACTCCGTCTCAAAAAAACAAAAACCGAAACAAGAACAAAAAACAAATTGAATGCCCTCAACTACCCTGTTTGGTGTAACTGCAGAGGTGGCCTGGAGATTAGTTCCAGGTGGCAGATGCCATGTGGCAGACAGACGCTGTTGTGGCAAGATACTGCCATGATGCCGGCGAAGCTGTGGTTTTTTTGGTGGCAGTGGGCACAGCAGGACCCCTTCCACACTGGGGAGCTGGGGCTGTGGGACCCCAGGTCTGGGCTCTGACGCTCTCAGCTGTCCAACCTCCCAGAGAGTC".replaceAll(" ", "");
    final String tmpl = "GAATTTTAAGTGGAATATTATCAGTGATGCATTTTTTACTGTAGTTTCTATCAATATAATTTTCACAATCGGCTGGGCGCAGTGGCTTACACCTGTAATACCAGCACTTTGGGAGGCTGAGGCGGGCGGATCACGGGGTCAGGGGATCAAAACCATCCTGGCTAACACGGTGAAACCCCGTCTCTACTAAAAATACAAAAAAAAATTAGCCGGGCATGGTGGCGGGCGCCTGTAGTCCCAGTTACTTGGGAGGCAGAGGCAAGAGAATGGTGTGAACCTGGGAGGCAGAGCTTGCAGTGAGCCAAGATCTCGCCACTGCACTCCAGCCTGGTCAACAGAGCGAGACTCCATCTCAAAAAAAAATAAATAAAAATAAATAAAATTTCACAATCAACTAAATTCTATAAAAATTCCATTGAGAGGCCGGGCGTGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGTGGATCATGAGGTCAGGAGATCGAGACCATCCTGGCTAACAAGGTGAAACCCCGTCTCTACTAAAAATACAAAAAATTAGCTGGGCGCGGTGGCGGGCACCTGTAGTCCCAGCTACTCGGGAGGCTGAGGCAGGAGAATGGCGTGAACCCGGGAAGCGGAGCTTGCAGTGAGCCGAGATTGCGCCACTGCAGTCTGCAGTCCGGCCTGGGCGACAGAGCGAGACTCCGTCTCAAAAAAAAAAAAAAAAAAAAAATTCCAT";
    final byte[] readbyte = DnaUtils.encodeString(read);
    final byte[] template = DnaUtils.encodeString(tmpl);
    final UnidirectionalEditDistance ed = getEditDistanceInstance(1, 1, 1, 0);
    final int[] actions = ed.calculateEditDistance(readbyte, readbyte.length, template, 30, 607, 250, true);
    if (actions != null) {
      assertEquals(read + " " + tmpl, 607, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
        assertEquals(30, actions[ActionsHelper.TEMPLATE_START_INDEX]);
    }
  }

  public void testAnotherNegSeedLen() {
    final String read = "    ATACTATATATTATATATATAACTATATTATATACATGTATATAATATATAACTATATATTATATATATAACTATTTATATACATGTATATAATATATAACTATATATTATATATATAACTATATTATATACATGTATATAATATATATAACTATATATTATATATATAACTATATTATATCATGTATATAATATATATAACTATATATTATATATAACTATATATTATATACATGTATATATTATATATAACTATATACTTATACATGTATATAATATATATAACTATTATATATAACTATATATTATATACATGTATATAACTATATTATATATAACTATATGTTATATACATGTATATAATATATAACTTATGTTATATATAACTATTCTATACGTGTATATAATATATAACTATATATTCTATACGTGTATATAATATATAACTATATATTCTATACGTGTATATAATATATAACTATATATTCTATACGTG".replaceAll(" ", "");
    final String tmpl = "TATATAACTATATATTATATATATAACTATATTATATACATGTATATAATATATAACTATATATTATATATATAACTATATTATATACATGTATATAATATATAACTATATATTATATATATAACTATATTATATACATGTATATAATATATAACTATATATTATATATATAACTATATTATATACATGTATATAATATATATAACTATATATTATATATATAACTATATTATATACATGTATATAATATATATAACTATATATTATATATAACTATATATTATATACATGTATATATTATATATAACTATATACTTTATACATGTATATAATATATATAACTATTATATATAACTATATATTATATACATGTATATAACTATATTATATATAACTATATGTTATATACATGTATATAATATATAACTATATGTTATATATAACTATTCTATACGTGTATATAATATATAACTATATATTCTATACGTGTATATAATATATAACTATATATTCTATACGTGTATAT";
    final byte[] readbyte = DnaUtils.encodeString(read);
    final byte[] template = DnaUtils.encodeString(tmpl);
    final UnidirectionalEditDistance ed = getEditDistanceInstance(1, 1, 1, 0);
    final int[] actions = ed.calculateEditDistance(readbyte, readbyte.length, template, 4, 50, MaxShiftUtils.calculateDefaultMaxShift(read.length()) * 2, true);
    if (actions != null) {
      fail();   //not enough coverage from seeds
    }
  }

  public void testYetAnotherNegSeedLen() {
    final String read = "ATCATCACCATCATCATCATCATCATCACCATCATCATCACCACCATCATCATCATCACCATCATCAACATCATCACCATCACCATCATCATCACCATCAC".replaceAll(" ", "");
    final String tmpl = "ATCATCACCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCACCATCATCATCATCATCATCATCATCAT";
    final byte[] readbyte = DnaUtils.encodeString(read);
    final byte[] template = DnaUtils.encodeString(tmpl);
    final UnidirectionalEditDistance ed = getEditDistanceInstance(1, 1, 1, 0);
    final int[] actions = ed.calculateEditDistance(readbyte, readbyte.length, template, 0, 50, MaxShiftUtils.calculateDefaultMaxShift(read.length()) * 2, true);
    if (actions != null) {
      fail();   //not enough coverage from seeds
    }
  }

  /*
   * this is a bad alignmend that seeded aligner creates (compare with gotoh)
   * public void testBadAlignment() {
    final String read = "    CACATACACACACACACACACACACATATATATATATATGTTTTCCCATGGTGCATTTAAAATAGAAACCAGTAATACAAAGAGAACAACAATACCAAAT".replaceAll(" ", "");;
    final String tmpl = "CACACACACACACACACACACACACACATATATATATATATATGTTTTCCCATGGTGCATTTAAAATAGAAACCAGTAATACAAAGAGAACAACAATACCAAATCCCCTGCATTTGGACATTT".replaceAll(" ", "");;
    final byte[] readbyte = DnaUtils.encodeString(read);
    final byte[] template = DnaUtils.encodeString(tmpl);
    final UnidirectionalEditDistance ed = getEditDistanceInstance(1, 1, false);
    final int[] actions = ed.calculateEditDistance(readbyte, readbyte.length, template, 4, 50, Utils.calculateDefaultMaxShift(read.length()) * 2, true);
    if (actions != null) {
      System.err.println(ActionsHelper.toString(actions));
    }
  }*/
}
