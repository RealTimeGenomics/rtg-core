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

import java.io.File;
import java.util.Random;

import org.apache.commons.lang.StringUtils;

import com.rtg.mode.DnaUtils;
import com.rtg.ngs.NgsParams;
import com.rtg.ngs.NgsParamsBuilder;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;


/**
 */
public class SingleIndelEditDistanceTest  extends AbstractUnidirectionalEditDistanceTest {

  private static final String[] DNA_STRINGS = new String[]{"A", "C", "G", "T"};

  private static String randomDNA(final long seed, final int length) {
    final Random rand = new Random(seed);
    final StringBuilder sb = new StringBuilder();
    for (int i = 0; i < length; i++) {
      final int dna = rand.nextInt(4);
      final String dnac = DNA_STRINGS[dna];
      sb.append(dnac);
    }
    return sb.toString();
  }
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
    final int m0 = SingleIndelEditDistance.m((byte) tb, (byte) rb, penalty, treatnsAsMismatches ? 9 : 0);
    assertEquals(exp, m0);
    final int m1 = SingleIndelEditDistance.m((byte) rb, (byte) tb, penalty, treatnsAsMismatches ? 9 : 0);
    assertEquals(exp, m1);
  }

  public void testDiagonal() throws Exception {
    final SingleIndelEditDistance ed = (SingleIndelEditDistance) getEditDistanceInstance(19, 1, 9, 0);
    final String exp = ""
        + "SingleIndelEditDistance" + " unknownsPenalty=0" + " substitution=9" + LS
        + " rLen=8" + " maxScore=2147483647" + " mMaxShift=9" + LS
        + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + LS
        + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + " 0" + LS
        + " firstMiss=0" + " lastMiss=7" + " diagScore=0" + " bestOffset=0" + LS
        ;
    final int[] actions = make("ACGTACGT", "ACGTACGT", 0, exp, ed);
    assertEquals(0, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
    assertEquals(0, actions[ActionsHelper.TEMPLATE_START_INDEX]);
    assertEquals("========", ActionsHelper.toString(actions));
  }

  public void testForwardPositive() throws Exception {
    final SingleIndelEditDistance ed = (SingleIndelEditDistance) getEditDistanceInstance(19, 1, 9, 0);
    for (int i = 1; i <= 9; i++) {
      checkForwardPositive(ed, i, null);
    }
  }

  public void testForwardPositive4() throws Exception {
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
    final String a = randomDNA(42, 20);
    final String b = randomDNA(101, 70);
    final String c = StringUtils.repeat("C", i);
    final String d = randomDNA(301, 30);
    final String e = randomDNA(21, 20);

    final String read = b + d;
    final String template = a + b + c + d + e;
    final int[] actions = make(read, template, a.length(), exp0, ed);
    final String exp = StringUtils.repeat("=", b.length()) + StringUtils.repeat("D", i) + StringUtils.repeat("=", d.length());
    assertEquals(exp, ActionsHelper.toString(actions));
    assertEquals(19 + i, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
    assertEquals(a.length(), actions[ActionsHelper.TEMPLATE_START_INDEX]);
  }

  public void testForwardNegative() throws Exception {
    final SingleIndelEditDistance ed = (SingleIndelEditDistance) getEditDistanceInstance(19, 1, 9, 0);
    for (int i = 1; i <= 9; i++) {
      checkForwardNegative(ed, i, null);
    }
  }

  public void testForwardNegative4() throws Exception {
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

  public void testForwardNegativeLong() throws Exception {
    final SingleIndelEditDistance ed = (SingleIndelEditDistance) getEditDistanceInstance(19, 1, 9, 0);
    for (int i = 1; i <= 9; i++) {
      checkForwardNegativeLong(ed, i, null);
    }
  }

  public void testForwardNegative4Long() throws Exception {
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
    final String a = randomDNA(42, 20);
    final String b = randomDNA(103, 70 - i);
    final String c = StringUtils.repeat("C", i);
    final String d = randomDNA(307, 30);
    final String e = randomDNA(71, 20);

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

  public void testReversePositive() throws Exception {
    final SingleIndelEditDistance ed = (SingleIndelEditDistance) getEditDistanceInstance(19, 1, 9, 4); //changed to true due to new off template allowance interfering with tests
    for (int i = 1; i <= 9; i++) {
      //final int i = 9;
      checkReversePositive(ed, i, null);
    }
  }

  public void testReversePositive4() throws Exception {
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

  public void testReverseNegative() throws Exception {
    final SingleIndelEditDistance ed = (SingleIndelEditDistance) getEditDistanceInstance(19, 1, 9, 0);
    for (int i = 1; i <= 9; i++) {
      checkReverseNegative(i, null, ed);
    }
  }

  public void testReverseNegative4() throws Exception {
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
    final String a = randomDNA(41, 20);
    final String b = randomDNA(503, 50);
    final String c = StringUtils.repeat("C", i);
    final String d = randomDNA(707, 50);
    final String e = randomDNA(61, 20);
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
    for (int i = 0; i < 10; i++) {
      getEditDistanceInstance(19, 1, 9, 0).calculateEditDistance(read, read.length, template, i, Integer.MAX_VALUE, 9, false);
    }
  }
  public void testThingInconsistent2() {
    final byte[] read = DnaUtils.encodeString("CCCAGGGAGCCGAGGGCGGTGGAGGTGGGGGCCCACCACCCAGCAGCAACTGTTGGACCCCGCAGCGCAGGAGTGGCAGGGGCCCGCTGCCGGGCAGGGGG");
    final byte[] template = DnaUtils.encodeString("GAAATAACCCCAGGGAGCCGAGGGCGGTGGAGGTGGGGGCCCACCACCCAGCAGCAACTGTTGGACCCCGCAGCGCAGGAGTGGCAGGGCCCGTCTGCAGTGCAGGAGTGGCAGGGCGTG");
    for (int i = 0; i < 18; i++) {
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

  public void testInsertAtStart() throws Exception {
    final byte[] read = DnaUtils.encodeString("      AAGACATTTTCTGATTTCATTATTCCAGTTGCCATTTCCTTAGTCAAGCCTTTCCTAATCTTCTTAACTAGGTCAAATTTACTATAATAGAGACTCATAGC".replaceAll(" ", ""));
    final byte[] template = DnaUtils.encodeString("GTATCCAGTAACCTGATTTCATTATTCCAGTTGCCATTTCCTTAGTCAAGCCTTTCCTAATCTTCTTAACTAGGTCAAATTTACTATAATAGAGACTCATAGCTTTTCTATCTCGAATTCATTGTCATGGTTATAATTTTACATATGTTT");
    runTestThing(12, read, template, "IIIIIIIIII===========================================================================================", 29);
    runTestThing(12 - 10, read, template, "IIIIIIIIII===========================================================================================", 29);
  }

  public void testInsertAtEnd() throws Exception {
    final byte[] read = DnaUtils.encodeString("      CTGATTTCATTATTCCAGTTGCCATTTCCTTAGTCAAGCCTTTCCTAATCTTCTTAACTAGGTCAAATTTACTATAATAGAGACTCATAGCAAGACAGGGG".replaceAll(" ", ""));
    final byte[] template = DnaUtils.encodeString("ACCTGATTTCATTATTCCAGTTGCCATTTCCTTAGTCAAGCCTTTCCTAATCTTCTTAACTAGGTCAAATTTACTATAATAGAGACTCATAGCTTTTCTATCTCGAATTC");
    runTestThing(2, read, template, "===========================================================================================IIIIIIIIII", 29);
  }

  public void testInsertAtEndWithMismatchAtStart() throws Exception {
    final byte[] read = DnaUtils.encodeString("      GTGATTTCATTATTCCAGTTGCCATTTCCTTAGTCAAGCCTTTCCTAATCTTCTTAACTAGGTCAAATTTACTATAATAGAGACTCATAGCAAGACAGGGG".replaceAll(" ", ""));
    final byte[] template = DnaUtils.encodeString("ACCTGATTTCATTATTCCAGTTGCCATTTCCTTAGTCAAGCCTTTCCTAATCTTCTTAACTAGGTCAAATTTACTATAATAGAGACTCATAGCTTTTCTATCTCGAATTC");
    runTestThing(2, read, template, "X==========================================================================================IIIIIIIIII", 38);
  }

  public void testNoMismatchAfterIndelAtEnd() throws Exception {
    final byte[] read = DnaUtils.encodeString("      TGGGCCTCTCGGTTCCCTCCCTATAATCTGTGGTTCACAGACAGTGGTGTGCACCAGAATCAGAACTGTGTGTTATGGCTGGGGCTTTTTTCCCCCTACAA".replaceAll(" ", ""));
    final byte[] template = DnaUtils.encodeString("CCTGGGCCTCTCGGTTCCCTCCCTATAATCTGTGGTTCACAGACAGTGGTGTGCACCAGAATCAGAACTGTGTGTTATGGCTGGTCCTTTTCGCCACTTACAAACAAGCC");
    runTestThing(2, read, template, "==================================================================================IIIIIIIIIIIIIIIIIII", 38);
  }


  public void testMatchBeforeInsert() throws Exception {
    final byte[] read = DnaUtils.encodeString("       CATGCATGCATACACACACACACACACACACACACACACAGAGTCTCTCTCTCTCTGTCTCTGTCTCTCCCTCTCTCACACACACACTCTCCATTCCCATC".replaceAll(" ", ""));
    final byte[] template = DnaUtils.encodeString("ACACACACACACACACACACACACACACACACACACACACACAGAGTCTCTCTCTCTCTGTCTCTGTCTCTCCCTCTCTCACACACACACTCTCCATTCCCATCATCCATGTCTGCTCCATCCTACAGGG");
    runTestThing(3, read, template, "=IIIIIIIIII==========================================================================================", 29);
  }
}
