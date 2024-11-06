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
package com.rtg.variant.bayes.multisample;

import com.rtg.mode.DNA;
import com.rtg.mode.DnaUtils;

import junit.framework.TestCase;

/**
 */
public class SingleNMerRepeatMeasurerTest extends TestCase {

  static byte[] makeRepeat(String repeat, int length) {
    final byte[] result = new byte[length];
    final byte[] rbytes = DnaUtils.encodeString(repeat);
    for (int i = 0; i < result.length; ++i) {
      result[i] = rbytes[i % rbytes.length];
    }
    return result;
  }

  public RepeatMeasurer getMeasurer(byte[] template) {
    return new SingleNMerRepeatMeasurer(template);
  }
  public RepeatMeasurer getMeasurer(byte[] template, int merLength) {
    return new SingleNMerRepeatMeasurer(template, merLength);
  }

  public void testMeasureRepeats() {
    //                 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18
    byte[] template = {2, 1, 3, 1, 2, 2, 2, 2, 1, 3, 2, 4, 2, 4, 2, 4, 2, 1, 3};
    RepeatMeasurer regionsA = getMeasurer(template);
    assertEquals(4, regionsA.measureRepeats(2, 9)); // Finds the 2  1-mer
    assertEquals(7, regionsA.measureRepeats(8, 18)); // Finds the 2,4  2-mer
    assertEquals(0, regionsA.measureRepeats(2, 18));

    //                     0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18
    template = new byte[] {2, 1, 3, 1, 2, 2, 1, 2, 2, 3, 2, 4, 2, 4, 2, 4, 3, 1, 3};
    regionsA =  getMeasurer(template);
    assertEquals(6, regionsA.measureRepeats(2, 9)); // Finds the 1,2,2  3-mer
    assertEquals(6, regionsA.measureRepeats(8, 18)); // Finds the 2,4  2-mer
    assertEquals(0, regionsA.measureRepeats(2, 18));

    //                     0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18
    template = new byte[] {2, 1, 3, 2, 1, 3, 2, 1, 3, 2, 1, 3, 2, 4, 2, 4, 3, 1, 3};
    regionsA =  getMeasurer(template);
    assertEquals(7, regionsA.measureRepeats(0, 7));
    assertEquals(8, regionsA.measureRepeats(0, 8));
    assertEquals(12, regionsA.measureRepeats(0, 12));
    //                     0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18
    template = new byte[] {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
    regionsA =  getMeasurer(template);
    assertEquals(18, regionsA.measureRepeats(0, 18));
  }

  public void testMeasureRepeatsJavadocExamples() {
    final byte[] templateA = DnaUtils.encodeString("accaccac");
    final byte[] templateB = DnaUtils.encodeString("accactacat");
    final byte[] templateC = DnaUtils.encodeString("accaaacac");
    final RepeatMeasurer regionsA =  getMeasurer(templateA);
    final RepeatMeasurer regionsB = getMeasurer(templateB);
    final RepeatMeasurer regionsC = getMeasurer(templateC);
//    * <code>accaccac</code> will count as a simple repeat of length 8 the third iteration is only partial but still counts.
    assertEquals(8, regionsA.measureRepeats(0, templateA.length + 1));
    assertEquals(8, regionsA.measureRepeats(0, templateA.length));

//    * <code>accac</code> is a repeat
    assertEquals(5, regionsA.measureRepeats(0, 5));

//    * <code>acca</code> is not a repeat unless followed by c
    assertEquals(4, regionsA.measureRepeats(0, 4)); // SRM3 doesn't look outside the range
    assertEquals(4, regionsB.measureRepeats(0, 4)); // SRM3 doesn't look outside the range
    assertEquals(2, regionsC.measureRepeats(0, 4));

//    * <code>aaa</code> is a repeat
    assertEquals(3, regionsC.measureRepeats(3, 6));

//    * <code>aca</code> is not a repeat (unless followed by c)
    assertEquals(0, regionsB.measureRepeats(6, 9));
    assertEquals(3, regionsC.measureRepeats(5, 8)); // SRM3 doesn't look outside the range

//    *<code>acc</code> is a repeat (cc)
    assertEquals(2, regionsC.measureRepeats(0, 3));
    assertEquals(3, regionsA.measureRepeats(0, 3));

    //    * <code>tt</code> is a repeat
    final byte[] template4 = DnaUtils.encodeString("tttcttc");
    final RepeatMeasurer regionsD = getMeasurer(template4);
    assertEquals(2, regionsD.measureRepeats(0, 2));
    assertEquals(2, regionsD.measureRepeats(4, 6));

    //ac is no repeat
    assertEquals(0, regionsC.measureRepeats(7, 9));
    assertEquals(2, regionsC.measureRepeats(5, 7));

    //    * <code>aatcagtt</code> is a repeat of aa, tt
    final byte[] template5 = DnaUtils.encodeString("aatcagtt");
    final RepeatMeasurer regionsE = getMeasurer(template5);
    assertEquals(0, regionsE.measureRepeats(0, template5.length)); // SRM3 does not chain together different repeat n-mers
  }

  public void testMeasureRealWorld() {
    final byte[] template = DNA.stringDNAtoByte("tttccaaaagagtttttcagagaaaaagagatgttgatacatttcaaaaatagga");
    assertEquals("agagaaaaagaga".length(), getMeasurer(template).measureRepeats(0, template.length - 1));
  }

  /*
  public void doExample(String template) {
    final byte[] tbytes = DNA.stringDNAtoByte(template);
    System.err.println(getMeasurer(tbytes,25).measureRepeats(0, tbytes.length) + "/" + tbytes.length);
  }

  // Tests a few examples
  public void testRealWorldExamples() {
    // Current is a little too aggro
    doExample("tagtgaggattggtgagtgtgaatgtgggtga"); // 1:1301941
    doExample("gcccgcccaccctggcttggc"); // 1:1887091
    doExample("cgtgggtgttaggttgtgggtaca"); // 1:1866605

    // Current is way too aggro
    doExample("acagaaacacacagacagaaacaaacacagagacacaga"); // 1:2258590
    doExample("gtctgtggtgtgtgtgtgcgtgtgtgtggtgtatgtgcatgtgtgtggtttgtgtggtgtgtgtagtgtgggtgtgtgtgtggt"); // 1:3326231
    doExample("cgtttcctataaatccttgggtggtgtgtgtgtatatat"); // 1:3366493
    doExample("ctgatgctgcttctcctcttcttcctcctcctcttcctcctctccttctcctcctcctctttttcctcctcctccttctccttcttcttcctttcttatttcttcttcctcttcctcttc"); // 1:4103800
    doExample("ttatatatattatatatatacaattataattatatattatatatacttataatatatataatatataat"); // 1:4160478
    doExample("ttccttctctccttccttctctcttcacctctctccttttctctctctttccttctctccctcctctta"); // 1:5306821
  }
  */

  public void test3Mer() {
    final byte[] template = makeRepeat("taa", 30);
    final RepeatMeasurer regionsA = getMeasurer(template, 3);
    assertEquals(template.length - 1, regionsA.measureRepeats(0, template.length - 1));
    assertEquals(template.length - 1, regionsA.measureRepeats(1, template.length));

    // Similarly...
    assertEquals(getMeasurer(makeRepeat("taa", 30)).measureRepeats(0, 30), getMeasurer(makeRepeat("ata", 30)).measureRepeats(0, 30));
    assertEquals(getMeasurer(makeRepeat("taa", 30)).measureRepeats(0, 30), getMeasurer(makeRepeat("atc", 30)).measureRepeats(0, 30));
    assertEquals(getMeasurer(makeRepeat("aaa", 30)).measureRepeats(0, 30), getMeasurer(makeRepeat("atc", 30)).measureRepeats(0, 30));
    assertEquals(getMeasurer(makeRepeat("taa", 30)).measureRepeats(0, 30), getMeasurer(makeRepeat("aat", 30)).measureRepeats(0, 30));
  }

  public void test4Mer() {
    byte[] template = makeRepeat("taaa", 30);
    RepeatMeasurer regionsA = getMeasurer(template, 4);
    assertEquals(template.length - 1, regionsA.measureRepeats(0, template.length - 1));
    assertEquals(template.length - 1, regionsA.measureRepeats(1, template.length));

    template = makeRepeat("ttaa", 30);
    regionsA = getMeasurer(template, 4);
    assertEquals(template.length - 1, regionsA.measureRepeats(0, template.length - 1));
    assertEquals(template.length - 1, regionsA.measureRepeats(1, template.length));

    template = makeRepeat("atta", 30);
    regionsA = getMeasurer(template, 4);
    assertEquals(template.length - 1, regionsA.measureRepeats(0, template.length - 1));
    assertEquals(template.length - 1, regionsA.measureRepeats(1, template.length));
  }

  public void test4Mer2() {
    assertEquals(getMeasurer(makeRepeat("atgc", 30), 4).measureRepeats(0, 30), getMeasurer(makeRepeat("taaa", 30), 4).measureRepeats(0, 30));
    assertEquals(getMeasurer(makeRepeat("atgc", 30), 4).measureRepeats(0, 30), getMeasurer(makeRepeat("ataa", 30), 4).measureRepeats(0, 30));
    assertEquals(getMeasurer(makeRepeat("atgc", 30), 4).measureRepeats(0, 30), getMeasurer(makeRepeat("aata", 30), 4).measureRepeats(0, 30));
    assertEquals(getMeasurer(makeRepeat("atgc", 30), 4).measureRepeats(0, 30), getMeasurer(makeRepeat("aaat", 30), 4).measureRepeats(0, 30));
    assertEquals(getMeasurer(makeRepeat("atgc", 30), 4).measureRepeats(0, 30), getMeasurer(makeRepeat("ttaa", 30), 4).measureRepeats(0, 30));
  }


}
