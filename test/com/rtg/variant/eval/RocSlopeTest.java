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

package com.rtg.variant.eval;

import static com.rtg.util.StringUtils.LS;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;

import junit.framework.TestCase;

/**
 */
public class RocSlopeTest extends TestCase {

  private static final String HEADER = "#posterior slope log-slope";

  private static final String ROC0 = (""
      + "10 1.0 1 seq1 1" + LS
      + "").replaceAll(" ", "\t");
  private static final String EXPECTED0 = (""
      + HEADER + LS
      + "10.00 1.00 0.000" + LS
      + "").replaceAll(" ", "\t");

  public void test0() throws IOException {
    final ByteArrayOutputStream baos = new ByteArrayOutputStream();
    try (PrintStream ps = new PrintStream(baos)) {
      RocSlope.writeSlope(new ByteArrayInputStream(ROC0.getBytes()), ps);
    }
    assertEquals(EXPECTED0, baos.toString());
  }

  private static final String ROC1 = ("#ignore me I'm a header" + LS
      + "10 1.0 0 seq1 1" + LS
      + "9 2.0 0 seq1 1" + LS
      + "9 3.0 0 seq1 1" + LS
      + "#ignore me too!" + LS
      + "9 3.0 1 seq1 1" + LS
      + "8 4.0 1 seq1 1" + LS
      + "7 5.0 1 seq1 1" + LS
      + "").replaceAll(" ", "\t");
  private static final String EXPECTED1 = (""
      + HEADER + LS
      + "9.00 2.00 0.301" + LS
      + "7.00 2.00 0.301" + LS
      + "").replaceAll(" ", "\t");

  public void test1() throws IOException {
    final ByteArrayOutputStream baos = new ByteArrayOutputStream();
    try (PrintStream ps = new PrintStream(baos)) {
      RocSlope.writeSlope(new ByteArrayInputStream(ROC1.getBytes()), ps);
    }
    assertEquals(EXPECTED1, baos.toString());
  }

  private static final String ROC2 = ("#ignore me I'm a header" + LS
      + "10 1.0 0 seq1 1" + LS
      + "9 2.0 0 seq1 1" + LS
      + "9 3.0 0 seq1 1" + LS
      + "#ignore me too!" + LS
      + "9 3.0 1 seq1 1" + LS
      + "8 4.0 1 seq1 1" + LS
      + "7 5.0 1 seq1 1" + LS
      + "6 5.0 2 seq1 1" + LS
      + "5 8.0 2 seq1 1" + LS
      + "4 8.0 4 seq1 1" + LS
      + "3 10.0 4 seq1 1" + LS
      + "2 10.0 5 seq1 1" + LS
      + "").replaceAll(" ", "\t");
  private static final String EXPECTED2 = (""
      + HEADER + LS
      + "9.00 2.50 0.398" + LS
      + "5.00 2.50 0.398" + LS
      + "5.00 1.00 0.000" + LS
      + "3.00 1.00 0.000" + LS
      + "").replaceAll(" ", "\t");

  public void test2() throws IOException {
    final ByteArrayOutputStream baos = new ByteArrayOutputStream();
    try (PrintStream ps = new PrintStream(baos)) {
      RocSlope.writeSlope(new ByteArrayInputStream(ROC2.getBytes()), ps);
    }
    assertEquals(EXPECTED2, baos.toString());
  }

  private static final String ROC3 = ("#ignore me I'm a header" + LS
      + "10 0.0 0 seq1 1" + LS
      + "9 0.0 1 seq1 1" + LS
      + "8 3.0 1 seq1 1" + LS
      + "7 3.0 2 seq1 1" + LS
      + "6 5.0 2 seq1 1" + LS
      + "5 5.0 3 seq1 1" + LS
      + "4 6.0 3 seq1 1" + LS
      + "").replaceAll(" ", "\t");
  private static final String EXPECTED3 = (""
      + HEADER + LS
      + "10.00 3.00 0.477" + LS
      + "8.00 3.00 0.477" + LS
      + "8.00 2.00 0.301" + LS
      + "6.00 2.00 0.301" + LS
      + "6.00 1.00 0.000" + LS
      + "4.00 1.00 0.000" + LS
      + "").replaceAll(" ", "\t");

  public void test3() throws IOException {
    final ByteArrayOutputStream baos = new ByteArrayOutputStream();
    try (PrintStream ps = new PrintStream(baos)) {
      RocSlope.writeSlope(new ByteArrayInputStream(ROC3.getBytes()), ps);
    }
    assertEquals(EXPECTED3, baos.toString());
  }

  private static final String ROC3_A = ("#ignore me I'm a header" + LS
      + "10 0.0 0 seq1 1" + LS
      + "8 3.0 1 seq1 1" + LS
      + "6 5.0 2 seq1 1" + LS
      + "4 6.0 3 seq1 1" + LS
      + "").replaceAll(" ", "\t");
  public void test3a() throws IOException {
    final ByteArrayOutputStream baos = new ByteArrayOutputStream();
    try (PrintStream ps = new PrintStream(baos)) {
      RocSlope.writeSlope(new ByteArrayInputStream(ROC3_A.getBytes()), ps);
    }
    assertEquals(EXPECTED3, baos.toString());
  }

}
