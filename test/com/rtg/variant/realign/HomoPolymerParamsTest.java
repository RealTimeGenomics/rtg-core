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

package com.rtg.variant.realign;

import static com.rtg.util.StringUtils.LS;

import java.io.IOException;
import java.io.Reader;
import java.io.StringReader;

import com.rtg.util.Utils;
import com.rtg.util.integrity.IntegralAbstract;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 */
public class HomoPolymerParamsTest extends TestCase {

  //typical case
  public void testCounts() throws IOException {
    final String inStr = ""
        + "A+T\t1\t1\t2\t3" + LS
        + "A+T\t3\t1\t2" + LS
        + "C+G\t3\t1\t2" + LS
        ;
    final Reader in = new StringReader(inStr);
    final int[][][] cnts = HomoPolymerParams.counts(in);
    final String exp = ""
        + "int[][][]" + LS
        + "[2]" + LS
        + "[0] int[][]" + LS
        + "  [4]" + LS
        + "  [0] null" + LS
        + "  [1] int[]" + LS
        + "    [3]" + LS
        + "    [0] 1, 2, 3" + LS
        + "  [2] null" + LS
        + "  [3] int[]" + LS
        + "    [2]" + LS
        + "    [0] 1, 2" + LS
        + "[1] int[][]" + LS
        + "  [4]" + LS
        + "  [0] null" + LS
        + "  [1] null" + LS
        + "  [2] null" + LS
        + "  [3] int[]" + LS
        + "    [2]" + LS
        + "    [0] 1, 2"
        ;
    final String str = IntegralAbstract.toString(cnts);
    //System.err.println(str);
    assertEquals(exp, str);
  }

  //syntax error in file - invalid tag
  public void testBadCounts1() {
    final String inStr = ""
        + "foo\t1\t1\t2\t3" + LS
        + "A+T\t3\t1\t2" + LS
        + "C+G\t3\t1\t2" + LS
        ;
    final Reader in = new StringReader(inStr);
    try {
      HomoPolymerParams.counts(in);
    } catch (final IOException e) {
      assertEquals("Invalid line in homopolymer calibration:foo 1 1 2 3", e.getMessage().replaceAll("\t", " "));
    }
  }

  //syntax error in file - too few fields
  public void testBadCounts2() {
    final String inStr = ""
        + "A+T" + LS
        + "A+T\t3\t1\t2" + LS
        + "C+G\t3\t1\t2" + LS
        ;
    final Reader in = new StringReader(inStr);
    try {
      HomoPolymerParams.counts(in);
    } catch (final IOException e) {
      assertEquals("Invalid line in homopolymer calibration:A+T", e.getMessage().replaceAll("\t", " "));
    }
  }

  public void testTranstions1a() {
    final double[] tr = HomoPolymerParams.transitions(SimplePossibility.SINGLETON, new int[] {}, false);
    assertEquals(0, tr.length);
  }

  public void testTranstions1b() {
    final double[] tr = HomoPolymerParams.transitions(SimplePossibility.SINGLETON, new int[] {0, 3, 2, 3, 0}, false);
    assertEquals(5, tr.length);
    assertEquals("[0.000, 0.375, 0.250, 0.375, 0.000]", Utils.realFormat(tr, 3));
  }

  //complement set true
  public void testTranstions1c() {
    final double[] tr = HomoPolymerParams.transitions(SimplePossibility.SINGLETON, new int[] {0, 3, 2, 3, 0}, true);
    assertEquals(5, tr.length);
    assertEquals("[1.000, 0.625, 0.750, 0.625, 1.000]", Utils.realFormat(tr, 3));
  }


  public void testTranstions3() {
    final double[][][] tr = HomoPolymerParams.transitions(SimplePossibility.SINGLETON, new int[][][] {{{1, 1}, {0, 1}}, {null, {1}}}, false);
    assertEquals(4, tr.length);
    final StringBuilder sb = new StringBuilder();
    for (int i = 0; i < 4; i++) {
      sb.append("nt:").append(i).append(LS);
      for (int t = 0; t < tr[i].length; t++) {
        sb.append("  t=").append(t).append(LS);
        sb.append("  ").append(Utils.realFormat(tr[i][t], 3)).append(LS);
      }
    }

    final String exp = ""
        + "nt:0" + LS
        + "  t=0" + LS
        + "  [0.500, 0.500]" + LS
        + "  t=1" + LS
        + "  [0.000, 1.000]" + LS
        + "nt:1" + LS
        + "  t=0" + LS
        + "  []" + LS
        + "  t=1" + LS
        + "  [1.000]" + LS
        + "nt:2" + LS
        + "  t=0" + LS
        + "  []" + LS
        + "  t=1" + LS
        + "  [1.000]" + LS
        + "nt:3" + LS
        + "  t=0" + LS
        + "  [0.500, 0.500]" + LS
        + "  t=1" + LS
        + "  [0.000, 1.000]" + LS
        ;
    final String str = sb.toString();
    //System.err.println(str);
    assertEquals(exp, str);
  }

  public void testTranstions3c() {
    final double[][][] tr = HomoPolymerParams.transitions(SimplePossibility.SINGLETON, new int[][][] {{{1, 1}, {0, 1}}, {null, {1}}}, true);
    assertEquals(4, tr.length);
    final StringBuilder sb = new StringBuilder();
    for (int i = 0; i < 4; i++) {
      sb.append("nt:").append(i).append(LS);
      for (int t = 0; t < tr[i].length; t++) {
        sb.append("  t=").append(t).append(LS);
        sb.append("  ").append(Utils.realFormat(tr[i][t], 3)).append(LS);
      }
    }

    final String exp = ""
        + "nt:0" + LS
        + "  t=0" + LS
        + "  [0.500, 0.500]" + LS
        + "  t=1" + LS
        + "  [1.000, 0.000]" + LS
        + "nt:1" + LS
        + "  t=0" + LS
        + "  []" + LS
        + "  t=1" + LS
        + "  [0.000]" + LS
        + "nt:2" + LS
        + "  t=0" + LS
        + "  []" + LS
        + "  t=1" + LS
        + "  [0.000]" + LS
        + "nt:3" + LS
        + "  t=0" + LS
        + "  [0.500, 0.500]" + LS
        + "  t=1" + LS
        + "  [1.000, 0.000]" + LS
        ;
    final String str = sb.toString();
    //System.err.println(str);
    assertEquals(exp, str);
  }

  public void test() throws IOException {
    final String inStr = ""
        + "A+T\t1\t1\t2\t3" + LS
        + "A+T\t3\t1\t2" + LS
        + "C+G\t1" + LS
        + "C+G\t3\t1\t2" + LS
        ;
    final Reader in = new StringReader(inStr);
    final HomoPolymerParams hpp = new HomoPolymerParams(SimplePossibility.SINGLETON, 2, 1, in);
    hpp.globalIntegrity();

    assertEquals(2, hpp.minReadRepeat());
    assertEquals(1, hpp.minTemplateRepeat());

    final StringBuilder sb = new StringBuilder();
    for (int i = 0; i < 4; i++) {
      sb.append("nt:").append(i).append(LS);
      for (int t = 1; t <= 3; t++) {
        sb.append("  t=").append(t).append(LS);
        for (int r = 1; r <= 3; r++) {
          sb.append("  ").append(Utils.realFormat(hpp.transition(i, t, r), 3)).append(":").append(Utils.realFormat(hpp.transitionC(i, t, r), 3));
        }
        sb.append(LS);
      }
    }
    final String exp = ""
        + "nt:0" + LS
        + "  t=1" + LS
        + "  0.000:1.000  0.333:0.667  0.500:0.500" + LS
        + "  t=2" + LS
        + "  0.000:1.000  0.000:1.000  0.000:1.000" + LS
        + "  t=3" + LS
        + "  0.000:1.000  0.667:0.333  0.000:1.000" + LS
        + "nt:1" + LS
        + "  t=1" + LS
        + "  0.000:1.000  0.000:1.000  0.000:1.000" + LS
        + "  t=2" + LS
        + "  0.000:1.000  0.000:1.000  0.000:1.000" + LS
        + "  t=3" + LS
        + "  0.000:1.000  0.667:0.333  0.000:1.000" + LS
        + "nt:2" + LS
        + "  t=1" + LS
        + "  0.000:1.000  0.000:1.000  0.000:1.000" + LS
        + "  t=2" + LS
        + "  0.000:1.000  0.000:1.000  0.000:1.000" + LS
        + "  t=3" + LS
        + "  0.000:1.000  0.667:0.333  0.000:1.000" + LS
        + "nt:3" + LS
        + "  t=1" + LS
        + "  0.000:1.000  0.333:0.667  0.500:0.500" + LS
        + "  t=2" + LS
        + "  0.000:1.000  0.000:1.000  0.000:1.000" + LS
        + "  t=3" + LS
        + "  0.000:1.000  0.667:0.333  0.000:1.000" + LS
        ;
    final String str = sb.toString();
    //System.err.println(str);
    assertEquals(exp, str);

  }
}
