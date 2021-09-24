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
package com.rtg.protein;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.Locale;

import com.rtg.AbstractTest;
import com.rtg.alignment.ActionsHelper;
import com.rtg.alignment.UnidirectionalEditDistance;
import com.rtg.mode.Protein;
import com.rtg.mode.ProteinScoringMatrix;
import com.rtg.mode.SequenceType;
import com.rtg.position.output.Blosum62;
import com.rtg.reader.CompressedMemorySequencesReader;
import com.rtg.util.InvalidParamsException;

/**
 */
public class GotohProteinEditDistanceTest extends AbstractTest {

  public UnidirectionalEditDistance getED(ProteinScoringMatrix matrix) {
    return new GotohProteinEditDistance(matrix);
  }

  private static String[] writeIgnore(final ProteinAlignmentResult res) throws IOException {
    final ByteArrayOutputStream bos = new ByteArrayOutputStream();
    try {
      res.write(bos);
    } finally {
      bos.close();
    }
    return bos.toString().split("\\t");
  }

  private CompressedMemorySequencesReader singleSequence(final byte[] seq) {
    return new CompressedMemorySequencesReader(new byte[][] {seq}, new String[] {"seq"}, new long[] {seq.length}, 0, 0, SequenceType.PROTEIN);
  }

  private void check(final String a, final String b, final String read, final String template, final String match, final int distance) throws InvalidParamsException, IOException {
    final byte[] s1 = Protein.encodeProteins(a.toUpperCase(Locale.getDefault()));
    final byte[] s2 = Protein.encodeProteins(b.toUpperCase(Locale.getDefault()));
    final UnidirectionalEditDistance f = getED(new ProteinScoringMatrix());
    final int[] r = f.calculateEditDistance(s1, s1.length, s2, 0, Integer.MAX_VALUE, 7, true);
    final SharedProteinResources resx = new SharedProteinResources(new ProteinScoringMatrix(), singleSequence(s2), singleSequence(s1), false);
    final ProteinAlignmentResult res = new ProteinAlignmentResult(resx, 0, 1, r, 0, true);
    final String[] p = writeIgnore(res);
    assertEquals(read, p[10]);
    assertEquals(template, p[9]);
    assertEquals(match, p[11]);
    assertEquals(distance, res.alignmentScore());
  }

  public void testProteinPerfectAlignment() throws InvalidParamsException, IOException {
    check("arndcqeghilkmfpstwyv", "arndcqeghilkmfpstwyv",
        "arndcqeghilkmfpstwyv",
        "arndcqeghilkmfpstwyv",
        "arndcqeghilkmfpstwyv", -116);
  }

  public void testProteinPerfectAlignmentWithUnknown() throws InvalidParamsException, IOException {
    check("arndcqeghilkmfpstwyv",
          "arndcqeghXlkmfpstwyv",
          "arndcqeghilkmfpstwyv",
          "arndcqeghxlkmfpstwyv",
          "arndcqegh lkmfpstwyv", -111);
  }

  public void testProteinSimilarAlignment() throws InvalidParamsException, IOException {
    check("arrnndqeehlkkkmmsstwyyvv",
          "sqkhsnrdqnirqeilansfhfim",
          "arrnndqeehlkkkmmsstwyyvv",
          "sqkhsnrdqnirqeilansfhfim",
          "++++++++++++++++++++++++", -35);
  }

  private int proteinExactMatchScore(final String s) {
    int score = 0;
    final Blosum62 matrix = Blosum62.SINGLETON;
    for (final byte b : Protein.encodeProteins(s.replace("+", "").replace(" ", "").replace("-", "").toUpperCase(Locale.getDefault()))) {
      score -= matrix.score(b, b);
    }
    return score;
  }

  public void testProteinSubAlignment() throws InvalidParamsException, IOException {
    final String r = "arndcqeghilkmfpstwyv";
    final String t = "arndcghilkmfpstwyv";
    final String u = "arndc--ghilkmfpstwyv";
    final String m = "arndc  ghilkmfpstwyv";
    check(r, t, r, u, m, proteinExactMatchScore(m) - (int) new ProteinScoringMatrix().getGap() - (2 * (int) new ProteinScoringMatrix().getExtend()));
  }

  public void testProteinindelAlignment() throws InvalidParamsException, IOException {
    final String r = "arndcghilkmfpstwyv";
    final String s = "arndc--ghilkmfpstwyv";
    final String t = "arndcqeghilkmfpstwyv";
    final String m = "arndc  ghilkmfpstwyv";
    check(r, t, s, t, m, proteinExactMatchScore(m) - (int) new ProteinScoringMatrix().getGap() - (2 * (int) new ProteinScoringMatrix().getExtend()));
  }

  public void testProteinSubAlignment2() throws InvalidParamsException, IOException {
    final String r = "arndcqeghilkmfpstwyv";
    final String t = "arndcewghilkmfpstwyv";
    final String m = "arndc+ ghilkmfpstwyv";
    check(r, t, r, t, m, proteinExactMatchScore(m) - 2 + 3);
  }

  public void testProteinSubAlignment3() throws InvalidParamsException, IOException {
    final String r = "arndcqeghilvywtspfmk";
    final String t = "arndcewghilkmfpstwyv";
    final String m = "arndc+ ghil  + s +  ";
    check(r, t, r, t, m, proteinExactMatchScore(m) - 2 + 3 + 2 + 1 - 1 + 1 + 1 - 1 + 1 + 2);
  }

  public void testProteinSubAlignment4() throws InvalidParamsException, IOException {
    final String r = "maqqrrggfkrrkkvdfiaa";
    final String m = "  +  r gf  r+k    + ";
    final String t =  "rrlvrrgfdwrrkdgtqsv*";
    // but Blast gives this (because it does not move the start pos?)
    //final String m = "    rrg   rrk       ";
    check(r, t, r, "x" + t.substring(0, t.length() - 1), m, -9);
  }

  public void testGoOffRightEndOfTemplate() throws InvalidParamsException, IOException {
    final UnidirectionalEditDistance ped = getED(new ProteinScoringMatrix());
    final byte[] read = {1, 2, 3};
    final byte[] template = {1, 2};
    final int[] alignment = ped.calculateEditDistance(read, read.length, template, 0, Integer.MAX_VALUE, 7, true);
    final SharedProteinResources res = new SharedProteinResources(new ProteinScoringMatrix(), singleSequence(template), singleSequence(read), false);
    final ProteinAlignmentResult r = new ProteinAlignmentResult(res, 0, 0, alignment, 0, true);
    final String[] p = writeIgnore(r);
    assertNotNull(r);
    assertEquals("*ax", p[9]);
    assertEquals("*ar", p[10]);
    assertEquals("*a ", p[11]);
    assertEquals(0, ActionsHelper.zeroBasedTemplateStart(alignment));
    assertEquals(-4, r.alignmentScore());
  }

  public void testGoOffLeftEndOfTemplate() throws InvalidParamsException, IOException {
    final UnidirectionalEditDistance ped = getED(new ProteinScoringMatrix());
    final byte[] read = {19, 20, 21};
    final byte[] template = {20, 21};
    final int[] alignment = ped.calculateEditDistance(read, read.length, template, -1, Integer.MAX_VALUE, 7, true);
    final SharedProteinResources res = new SharedProteinResources(new ProteinScoringMatrix(), singleSequence(template), singleSequence(read), false);
    final ProteinAlignmentResult r = new ProteinAlignmentResult(res, 0, 0, alignment, 0, true);
    final String[] p = writeIgnore(r);
    assertEquals("xyv", p[9]);
    assertEquals("wyv", p[10]);
    assertEquals(" yv", p[11]);
    assertEquals(-1, ActionsHelper.zeroBasedTemplateStart(alignment));
    assertEquals(-9, r.alignmentScore());
  }
}
