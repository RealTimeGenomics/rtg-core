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

import static com.rtg.util.StringUtils.LS;

import java.io.File;

import org.apache.commons.lang.StringUtils;

import com.rtg.mode.DnaUtils;
import com.rtg.ngs.NgsParams;
import com.rtg.ngs.NgsParamsBuilder;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.RandomDna;


/**
 */
public class SingleIndelEditDistanceTest  extends AbstractUnidirectionalEditDistanceTest {

  @Override
  protected UnidirectionalEditDistance getEditDistanceInstance(int gapOpenPenalty, int gapExtendPenalty, int substitutionPenalty, int unknownsPenalty) {
        final NgsParams params = new NgsParamsBuilder().gapOpenPenalty(gapOpenPenalty).gapExtendPenalty(gapExtendPenalty).substitutionPenalty(substitutionPenalty).unknownsPenalty(unknownsPenalty).singleIndelPenalties(null).create();
        return new SingleIndelEditDistance(params, 1000);
  }

  private SingleIndelEditDistance getTableEDInstance(String table, int maxReadLen) throws Exception {
    try (TestDirectory dir = new TestDirectory("siedt")) {
      final File tableFile = FileUtils.stringToFile(table, new File(dir, "ttt.properties"));
      final NgsParams params = new NgsParamsBuilder().singleIndelPenalties(tableFile.getPath()).create();
      return new SingleIndelEditDistance(params, maxReadLen);
    }
  }

  public void testTableLoading2() throws Exception {
    String table = ""
                           + "error_ins_penalty_extension_slope=1.0\n"
                           + "error_del_penalty_extension_slope=0.5\n"
                           + "error_ins_penalty=5,6,6,6,6,6\n" //...        7,   8,   9   10
                           + "error_del_penalty=4,6,6,6,6\n"  //... 6.5, 7.0, 7.5, 8.0, 8.5
            ;
    try {
      getTableEDInstance(table, 15);
      fail();
    } catch (NoTalkbackSlimException ntse) {
      TestUtils.containsAll(ntse.getMessage(), "Could not parse substitution penalty in");
    }

    table += "error_snp_penalty=9\n";

    try {
      getTableEDInstance(table, 15);
      fail();
    } catch (NoTalkbackSlimException ntse) {
      TestUtils.containsAll(ntse.getMessage(), "Could not parse unknowns penalty in");
    }

    table += "error_unknowns_penalty=5\n";
    final SingleIndelEditDistance ed = getTableEDInstance(table, 10);

    TestUtils.assertEquals(new int[]{4,  5,  6, 6,  6, 6,  6, 6,  6, 6,  6, /* now extrapolation: */ 6,  7,  7,  7, 8,  8,  8, 9, 10}, ed.getIndelPenalties());
    TestUtils.assertEquals(new int[]{-1, 1, -2, 2, -3, 3, -4, 4, -5, 5, -6,                          6, -7,  7, -8, 8, -9, -10, 9, 10}, ed.getIndelOffsets());
  }

  public void testM() {
    checkm(1, 1, true, 0);
    checkm(1, 1, false, 0);
    checkm(2, 1, true, 9);
    checkm(2, 1, false, 9);
    checkm(4, 4, true, 0);
    checkm(4, 4, false, 0);
    checkm(1, 0, true, 9);
    checkm(1, 0, false, 0);
    checkm(0, 0, true, 9);
    checkm(0, 0, false, 0);
  }

  private void checkm(final int rb, final int tb, final boolean treatnsAsMismatches, final int exp) {
    final int penalty = 9;
    final int m0 = SingleIndelEditDistance.m((byte) tb, (byte) rb, 0, penalty, treatnsAsMismatches ? 9 : 0);
    assertEquals(exp, m0);
    final int m1 = SingleIndelEditDistance.m((byte) rb, (byte) tb, 0, penalty, treatnsAsMismatches ? 9 : 0);
    assertEquals(exp, m1);
  }

  public void testDiagonal() {
    final SingleIndelEditDistance ed = (SingleIndelEditDistance) getEditDistanceInstance(19, 1, 9, 0);
    final String exp = ""
        + "SingleIndelEditDistance" + " unknownsPenalty=0" + " substitution=9" + LS
        + " rLen=8" + " maxScore=2147483647" + " mMaxShift=9" + LS
        + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + LS
        + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + LS
        + " firstMiss=-1" + " lastMiss=-1" + " diagScore=0" + " bestOffset=0" + LS
        ;
    final int[] actions = make("ACGTACGT", "ACGTACGT", 0, exp, ed);
    assertEquals(0, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
    assertEquals(0, actions[ActionsHelper.TEMPLATE_START_INDEX]);
    assertEquals("========", ActionsHelper.toString(actions));
  }

  public void testForwardPositive() {
    final SingleIndelEditDistance ed = (SingleIndelEditDistance) getEditDistanceInstance(19, 1, 9, 0);
    for (int i = 1; i <= 9; ++i) {
      checkForwardPositive(ed, i, null);
    }
  }

  public void testForwardPositive4() {
    final SingleIndelEditDistance ed = (SingleIndelEditDistance) getEditDistanceInstance(19, 1, 9, 0);
    final String exp = ""
        + "SingleIndelEditDistance" + " unknownsPenalty=0" + " substitution=9" + LS
        + " rLen=100" + " maxScore=2147483647" + " mMaxShift=9" + LS
        + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 9" + " 9" + " 0" + " 9" + " 9" + " 9" + " 0" + " 9" + " 0" + " 9" + " 0" + " 9" + " 0" + " 9" + " 0" + " 0" + " 9" + " 9" + " 9" + " 9" + " 9" + " 9" + " 0" + " 9" + " 9" + " 9" + " 9" + " 9" + " 9" + " 9" + LS
        + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 9" + " 18" + " 18" + " 27" + " 36" + " 45" + " 45" + " 54" + " 54" + " 63" + " 63" + " 72" + " 72" + " 81" + " 81" + " 81" + " 90" + " 99" + " 108" + " 117" + " 126" + " 135" + " 135" + " 144" + " 153" + " 162" + " 171" + " 180" + " 189" + LS
        + " firstMiss=70" + " lastMiss=99" + " diagScore=198" + " bestOffset=4" + LS
        + " bestPosn=70" + " bestScore=23" + " bestForward=true" + LS
        ;
    checkForwardPositive(ed, 4, exp);
  }

  private void checkForwardPositive(final SingleIndelEditDistance ed, int i, final String exp0) {
    //    final String tins = "AAAA" + StringUtils.repeat("C", i) + "GTAT";
    //    final int[] actions = make("AAAAGTAT", tins, 0, exp0, ed);
    //    final String exp = "====" + StringUtils.repeat("D", i) + "====";
    final String a = RandomDna.random(20, 42);
    final String b = RandomDna.random(70, 101);
    final String c = StringUtils.repeat("C", i);
    final String d = RandomDna.random(30, 301);
    final String e = RandomDna.random(20, 21);

    final String read = b + d;
    final String template = a + b + c + d + e;
    final int[] actions = make(read, template, a.length(), exp0, ed);
    final String exp = StringUtils.repeat("=", b.length()) + StringUtils.repeat("D", i) + StringUtils.repeat("=", d.length());
    assertEquals(exp, ActionsHelper.toString(actions));
    assertEquals(19 + i, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
    assertEquals(a.length(), actions[ActionsHelper.TEMPLATE_START_INDEX]);
  }

  public void testForwardNegative() {
    final SingleIndelEditDistance ed = (SingleIndelEditDistance) getEditDistanceInstance(19, 1, 9, 0);
    for (int i = 1; i <= 9; ++i) {
      checkForwardNegative(ed, i, null);
    }
  }

  public void testForwardNegative4() {
    final SingleIndelEditDistance ed = (SingleIndelEditDistance) getEditDistanceInstance(19, 1, 9, 0);
    final String exp = ""
        + "SingleIndelEditDistance" + " unknownsPenalty=0" + " substitution=9" + LS
        + " rLen=12" + " maxScore=2147483647" + " mMaxShift=9" + LS
        + " 0" + " 0" + " 0" + " 0" + " 9" + " 9" + " 9" + " 9" + " 0" + " 0" + " 0" + " 0" + LS
        + " 0" + " 0" + " 0" + " 0" + " 0" + " 9" + " 18" + " 27" + " 36" + " 36" + " 36" + " 36" + LS
        + " firstMiss=4" + " lastMiss=7" + " diagScore=36" + " bestOffset=-4" + LS
        + " bestPosn=8" + " bestScore=23" + " bestForward=true" + LS
        ;
    checkForwardNegative(ed, 4, exp);
  }

  private void checkForwardNegative(final SingleIndelEditDistance ed, int i, final String exp0) {
    final String read = "AAAA" + StringUtils.repeat("C", i) + "GTAT";
    final int[] actions = make(read, "AAAA" + StringUtils.repeat("GTAT", 4), 0, exp0, ed);
    final String exp = "====" + StringUtils.repeat("I", i) + "====";
    assertEquals(exp, ActionsHelper.toString(actions));
    assertEquals(19 + i, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
    assertEquals(0, actions[ActionsHelper.TEMPLATE_START_INDEX]);
  }

  public void testForwardNegativeLong() {
    final SingleIndelEditDistance ed = (SingleIndelEditDistance) getEditDistanceInstance(19, 1, 9, 0);
    for (int i = 1; i <= 9; ++i) {
      checkForwardNegativeLong(ed, i, null);
    }
  }

  public void testForwardNegative4Long() {
    final SingleIndelEditDistance ed = (SingleIndelEditDistance) getEditDistanceInstance(19, 1, 9, 0);
    final String exp = ""
        + "SingleIndelEditDistance" + " unknownsPenalty=0" + " substitution=9" + LS
        + " rLen=100" + " maxScore=2147483647" + " mMaxShift=9" + LS
        + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 9" + " 9" + " 0" + " 9" + " 9" + " 9" + " 9" + " 0" + " 9" + " 9" + " 0" + " 9" + " 9" + " 0" + " 9" + " 0" + " 9" + " 9" + " 0" + " 9" + " 9" + " 9" + " 9" + " 9" + " 9" + " 9" + " 9" + " 9" + " 9" + " 9" + " 0" + " 0" + " 9" + " 0" + LS
        + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 9" + " 18" + " 18" + " 27" + " 36" + " 45" + " 54" + " 54" + " 63" + " 72" + " 72" + " 81" + " 90" + " 90" + " 99" + " 99" + " 108" + " 117" + " 117" + " 126" + " 135" + " 144" + " 153" + " 162" + " 171" + " 180" + " 189" + " 198" + " 207" + " 216" + " 216" + " 216" + " 225" + LS
        + " firstMiss=66" + " lastMiss=98" + " diagScore=225" + " bestOffset=-4" + LS
        + " bestPosn=70" + " bestScore=23" + " bestForward=true" + LS
        ;
    checkForwardNegativeLong(ed, 4, exp);
  }

  private void checkForwardNegativeLong(final SingleIndelEditDistance ed, int i, final String exp0) {
    final String a = RandomDna.random(20, 42);
    final String b = RandomDna.random(70 - i, 103);
    final String c = StringUtils.repeat("C", i);
    final String d = RandomDna.random(30, 307);
    final String e = RandomDna.random(20, 71);

    final String read = b + c + d;
    final String template = a + b + d + e;
    final int[] actions = make(read, template, a.length(), exp0, ed);
    final String exp = StringUtils.repeat("=", b.length()) + StringUtils.repeat("I", i) + StringUtils.repeat("=", d.length());
    assertEquals(exp, ActionsHelper.toString(actions));
    assertEquals(19 + i, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
    assertEquals(a.length(), actions[ActionsHelper.TEMPLATE_START_INDEX]);
  }

  private static final String TRICKY = "GTAGGTTATAGATTAGAT";
  private static String sub(final int i) {
    return TRICKY.substring(0, i);
  }

  public void testReversePositive() {
    final SingleIndelEditDistance ed = (SingleIndelEditDistance) getEditDistanceInstance(19, 1, 9, 4); //changed to true due to new off template allowance interfering with tests
    for (int i = 1; i <= 9; ++i) {
      //final int i = 9;
      checkReversePositive(ed, i, null);
    }
  }

  public void testReversePositive4() {
    final SingleIndelEditDistance ed = (SingleIndelEditDistance) getEditDistanceInstance(19, 1, 9, 0);
    final String exp0 = ""
        + "SingleIndelEditDistance" + " unknownsPenalty=0" + " substitution=9" + LS
        + " rLen=25" + " maxScore=2147483647" + " mMaxShift=9" + LS
        + " 9" + " 9" + " 9" + " 9" + " 0" + " 0" + " 9" + " 9" + " 9" + " 9" + " 9" + " 9" + " 9" + " 9" + " 9" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + LS
        + " 0" + " 9" + " 18" + " 27" + " 36" + " 36" + " 36" + " 45" + " 54" + " 63" + " 72" + " 81" + " 90" + " 99" + " 108" + " 117" + " 117" + " 117" + " 117" + " 117" + " 117" + " 117" + " 117" + " 117" + " 117" + LS
        + " firstMiss=0" + " lastMiss=14" + " diagScore=117" + " bestOffset=4" + LS
        + " bestPosn=10" + " bestScore=23" + " bestForward=false" + LS
        ;
    checkReversePositive(ed, 4, exp0);
  }

  private void checkReversePositive(final SingleIndelEditDistance ed, int i, final String exp0) {

    final String a = StringUtils.repeat("C", i);
    final String b = sub(15 - i);
    final String c = StringUtils.repeat("C", i);
    final String d = "GTATGTATAA";
    final String template = a + b + d;
    final String read = b + c + d;
    final int[] actions = make(read, template, 0, exp0, ed);
    final String exp = StringUtils.repeat("=", b.length()) + StringUtils.repeat("I", c.length()) + StringUtils.repeat("=", d.length());
    assertEquals(exp, ActionsHelper.toString(actions));
    assertEquals(19 + i, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
    assertEquals(i, actions[ActionsHelper.TEMPLATE_START_INDEX]);
  }

  public void testReverseNegative() {
    final SingleIndelEditDistance ed = (SingleIndelEditDistance) getEditDistanceInstance(19, 1, 9, 0);
    for (int i = 1; i <= 9; ++i) {
      checkReverseNegative(i, null, ed);
    }
  }

  public void testReverseNegative4() {
    final SingleIndelEditDistance ed = (SingleIndelEditDistance) getEditDistanceInstance(19, 1, 9, 0);
    final String exp0 = ""
        + "SingleIndelEditDistance" + " unknownsPenalty=0" + " substitution=9" + LS
        + " rLen=100" + " maxScore=2147483647" + " mMaxShift=9" + LS
        + " 0" + " 9" + " 0" + " 9" + " 9" + " 9" + " 9" + " 9" + " 9" + " 0" + " 9" + " 9" + " 0" + " 9" + " 9" + " 9" + " 0" + " 9" + " 9" + " 9" + " 9" + " 9" + " 9" + " 0" + " 9" + " 9" + " 9" + " 9" + " 9" + " 9" + " 9" + " 9" + " 9" + " 0" + " 9" + " 9" + " 0" + " 9" + " 9" + " 9" + " 9" + " 9" + " 0" + " 9" + " 0" + " 0" + " 9" + " 9" + " 9" + " 9" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + LS
        + " 0" + " 0" + " 9" + " 9" + " 18" + " 27" + " 36" + " 45" + " 54" + " 63" + " 63" + " 72" + " 81" + " 81" + " 90" + " 99" + " 108" + " 108" + " 117" + " 126" + " 135" + " 144" + " 153" + " 162" + " 162" + " 171" + " 180" + " 189" + " 198" + " 207" + " 216" + " 225" + " 234" + " 243" + " 243" + " 252" + " 261" + " 261" + " 270" + " 279" + " 288" + " 297" + " 306" + " 306" + " 315" + " 315" + " 315" + " 324" + " 333" + " 342" + " 351" + " 351" + " 351" + " 351" + " 351" + " 351" + " 351" + " 351" + " 351" + " 351" + " 351" + " 351" + " 351" + " 351" + " 351" + " 351" + " 351" + " 351" + " 351" + " 351" + " 351" + " 351" + " 351" + " 351" + " 351" + " 351" + " 351" + " 351" + " 351" + " 351" + " 351" + " 351" + " 351" + " 351" + " 351" + " 351" + " 351" + " 351" + " 351" + " 351" + " 351" + " 351" + " 351" + " 351" + " 351" + " 351" + " 351" + " 351" + " 351" + " 351" + LS
        + " firstMiss=1" + " lastMiss=49" + " diagScore=351" + " bestOffset=-4" + LS
        + " bestPosn=49" + " bestScore=23" + " bestForward=false" + LS
        ;
    checkReverseNegative(4, exp0, ed);
  }

  private void checkReverseNegative(int i, final String exp0, final SingleIndelEditDistance ed) {
    final String a = RandomDna.random(20, 41);
    final String b = RandomDna.random(50, 503);
    final String c = StringUtils.repeat("C", i);
    final String d = RandomDna.random(50, 707);
    final String e = RandomDna.random(20, 61);
    //    final String a = "";
    //    final String b = sub(15 - i);
    //    final String c = StringUtils.repeat("C", i);
    //    final String d = "GTATGTAT";
    //    final String e = "";
    final String template = a + b + c + d + e;
    final String read = b + d;
    final int[] actions = make(read, template, a.length() + i, exp0, ed);
    final String exp = StringUtils.repeat("=", b.length()) + StringUtils.repeat("D", c.length()) + StringUtils.repeat("=", d.length());
    assertEquals(exp, ActionsHelper.toString(actions));
    assertEquals(19 + i, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
    assertEquals(a.length(), actions[ActionsHelper.TEMPLATE_START_INDEX]);
  }

  private int[] make(final String readStr, final String templateStr, final int tStart, final String exp, final SingleIndelEditDistance ed) {
    //System.err.println(" read       =" + readStr);
    //System.err.println(" templateStr=" + templateStr);
    final byte[] read = DnaUtils.encodeString(readStr);
    final byte[] template = DnaUtils.encodeString(templateStr);

    final int[] actions = ed.calculateEditDistance(read, read.length, template, tStart, Integer.MAX_VALUE, 9, false);
    //System.err.println(ActionsHelper.toString(actions));
    ed.globalIntegrity();
    if (exp != null) {
      assertEquals(exp, ed.toString());
    }
    return actions;
  }

  public void testIndelOffTemplateEnd() {
    final String readStr =     "AAATCCACATTATAT     AAATCTAACTATATCATCATACAATACATACATAG".replaceAll(" ", "");
    final byte[] read = DnaUtils.encodeString(readStr);
    final String templateStr = "AAATCCACATTATATGGGGGAAATCTAACTATATCATCATACAATACATACATA";
    final byte[] template = DnaUtils.encodeString(templateStr);

    final int[] actions = getEditDistanceInstance(19, 1, 9, 0).calculateEditDistance(read, read.length, template, 0, Integer.MAX_VALUE, 9, false);
    final String exp = "===============DDDDD==================================X"; //regression
    assertEquals(exp, ActionsHelper.toString(actions));
    assertEquals(0, actions[ActionsHelper.TEMPLATE_START_INDEX]);
    //System.err.println(ActionsHelper.toString(actions)); //not sure what the intended behaviour is here, return null? "best alignment" that won't leave the template?
  }

  //these two tests (testInconsistent) trigger an assert that fails before a recent bugfix (Tues 1 Oct 2013, a9380b0)
  public void testThingInconsistent() {
    final byte[] read = DnaUtils.encodeString("GGGTGGGTCGGGTGGGGTGTGGGGGGGGCTGGTGAGGTGTGGGGGGGGCTGGTCAGGCGTGGGCTGGGCTGGTCAGGCGTGGGGTCAGCCTGGCCCCTGCC");
    final byte[] template = DnaUtils.encodeString("CATGGGGCCGGCTGGTCAGGCGTGGGGCGGGCTGGTCAGGCGTGGGGCGGGCTGGTCAGGCGTGGGCTGGGCTGGTCAGGCGTGGGGTCAGCCTGGCCCCTGCCCGCTGG");
    for (int i = 0; i < 10; ++i) {
      getEditDistanceInstance(19, 1, 9, 0).calculateEditDistance(read, read.length, template, i, Integer.MAX_VALUE, 9, false);
    }
  }
  public void testThingInconsistent2() {
    final byte[] read = DnaUtils.encodeString("CCCAGGGAGCCGAGGGCGGTGGAGGTGGGGGCCCACCACCCAGCAGCAACTGTTGGACCCCGCAGCGCAGGAGTGGCAGGGGCCCGCTGCCGGGCAGGGGG");
    final byte[] template = DnaUtils.encodeString("GAAATAACCCCAGGGAGCCGAGGGCGGTGGAGGTGGGGGCCCACCACCCAGCAGCAACTGTTGGACCCCGCAGCGCAGGAGTGGCAGGGCCCGTCTGCAGTGCAGGAGTGGCAGGGCGTG");
    for (int i = 0; i < 18; ++i) {
      getEditDistanceInstance(19, 1, 9, 0).calculateEditDistance(read, read.length, template, i, Integer.MAX_VALUE, 9, false);
    }
  }


  private void runTestThing(int position, byte[] read, byte[] tmpl, String actionString, int alignmentScore) {
    final UnidirectionalEditDistance ed = getEditDistanceInstance(19, 1, 9, 0);
    final int[] actions = ed.calculateEditDistance(read, read.length, tmpl, position, 255, 50, true);
    if (actions == null) {
//        System.err.println("Null");
      assertNull(actionString);
    } else {
//        System.err.println(ActionsHelper.toString(actions));
      assertNotNull(actionString);
      assertEquals(actionString, ActionsHelper.toString(actions));
      assertEquals(alignmentScore, ActionsHelper.alignmentScore(actions));
    }
  }

  public void testInsertAtStart() {
    final byte[] read = DnaUtils.encodeString("      AAGACATTTTCTGATTTCATTATTCCAGTTGCCATTTCCTTAGTCAAGCCTTTCCTAATCTTCTTAACTAGGTCAAATTTACTATAATAGAGACTCATAGC".replaceAll(" ", ""));
    final byte[] template = DnaUtils.encodeString("GTATCCAGTAACCTGATTTCATTATTCCAGTTGCCATTTCCTTAGTCAAGCCTTTCCTAATCTTCTTAACTAGGTCAAATTTACTATAATAGAGACTCATAGCTTTTCTATCTCGAATTCATTGTCATGGTTATAATTTTACATATGTTT");
    runTestThing(12, read, template, "IIIIIIIIII===========================================================================================", 29);
    runTestThing(12 - 10, read, template, "IIIIIIIIII===========================================================================================", 29);
  }

  public void testInsertAtEnd() {
    final byte[] read = DnaUtils.encodeString("      CTGATTTCATTATTCCAGTTGCCATTTCCTTAGTCAAGCCTTTCCTAATCTTCTTAACTAGGTCAAATTTACTATAATAGAGACTCATAGCAAGACAGGGG".replaceAll(" ", ""));
    final byte[] template = DnaUtils.encodeString("ACCTGATTTCATTATTCCAGTTGCCATTTCCTTAGTCAAGCCTTTCCTAATCTTCTTAACTAGGTCAAATTTACTATAATAGAGACTCATAGCTTTTCTATCTCGAATTC");
    runTestThing(2, read, template, "===========================================================================================IIIIIIIIII", 29);
  }

  public void testInsertAtEndWithMismatchAtStart() {
    final byte[] read = DnaUtils.encodeString("      GTGATTTCATTATTCCAGTTGCCATTTCCTTAGTCAAGCCTTTCCTAATCTTCTTAACTAGGTCAAATTTACTATAATAGAGACTCATAGCAAGACAGGGG".replaceAll(" ", ""));
    final byte[] template = DnaUtils.encodeString("ACCTGATTTCATTATTCCAGTTGCCATTTCCTTAGTCAAGCCTTTCCTAATCTTCTTAACTAGGTCAAATTTACTATAATAGAGACTCATAGCTTTTCTATCTCGAATTC");
    runTestThing(2, read, template, "X==========================================================================================IIIIIIIIII", 38);
  }

  public void testNoMismatchAfterIndelAtEnd() {
    final byte[] read = DnaUtils.encodeString("      TGGGCCTCTCGGTTCCCTCCCTATAATCTGTGGTTCACAGACAGTGGTGTGCACCAGAATCAGAACTGTGTGTTATGGCTGGGGCTTTTTTCCCCCTACAA".replaceAll(" ", ""));
    final byte[] template = DnaUtils.encodeString("CCTGGGCCTCTCGGTTCCCTCCCTATAATCTGTGGTTCACAGACAGTGGTGTGCACCAGAATCAGAACTGTGTGTTATGGCTGGTCCTTTTCGCCACTTACAAACAAGCC");
    runTestThing(2, read, template, "==================================================================================IIIIIIIIIIIIIIIIIII", 38);
  }


  public void testMatchBeforeInsert() {
    final byte[] read = DnaUtils.encodeString("       CATGCATGCATACACACACACACACACACACACACACACAGAGTCTCTCTCTCTCTGTCTCTGTCTCTCCCTCTCTCACACACACACTCTCCATTCCCATC".replaceAll(" ", ""));
    final byte[] template = DnaUtils.encodeString("ACACACACACACACACACACACACACACACACACACACACACAGAGTCTCTCTCTCTCTGTCTCTGTCTCTCCCTCTCTCACACACACACTCTCCATTCCCATCATCCATGTCTGCTCCATCCTACAGGG");
    runTestThing(3, read, template, "=IIIIIIIIII==========================================================================================", 29);
  }
}
