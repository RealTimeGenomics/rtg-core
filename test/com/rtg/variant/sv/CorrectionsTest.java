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

package com.rtg.variant.sv;

import static com.rtg.util.StringUtils.LS;

import java.io.ByteArrayInputStream;
import java.io.IOException;

import junit.framework.TestCase;

/**
 */
public class CorrectionsTest extends TestCase {

  public void test() throws IOException {
    final String coStr = ""
      + "# a comment" + LS
      + LS
      + "seq1 0 1.0" + LS
      + "seq1 1 2.0" + LS
      + "seq1 2 0.25" + LS
      + "seq1 5 1.5" + LS
      ;
    final Corrections co = new Corrections(new ByteArrayInputStream(coStr.replace(' ', '\t').getBytes()), 1);
    assertEquals(9, co.length());
    assertEquals(1.0, co.correction(0));
    assertEquals(2.0, co.correction(1));
    assertEquals(0.25, co.correction(2));
    assertEquals(0.0, co.correction(3));
    assertEquals(0.0, co.correction(4));
    assertEquals(1.5, co.correction(5));
  }

  public void testBad1() {
    final String coStr = ""
      + " # a comment" + LS
      + "seq1 0 1.0" + LS
      + "seq1 1 2.0" + LS
      + "seq1 2 0.25" + LS
      + "seq1 5 1.5" + LS
      ;
    try {
      new Corrections(new ByteArrayInputStream(coStr.replace(' ', '\t').getBytes()), 1);
      fail();
    } catch (final Exception e) {
      assertEquals(" # a comment", e.getMessage().replace('\t', ' '));
    }
  }

  public void testBad2() {
    final String coStr = ""
      + "# a comment" + LS
      + "seq1 0" + LS
      + "seq1 1 2.0" + LS
      + "seq1 2 0.25" + LS
      + "seq1 5 1.5" + LS
      ;
    try {
      new Corrections(new ByteArrayInputStream(coStr.replace(' ', '\t').getBytes()), 1);
      fail();
    } catch (final Exception e) {
      assertEquals("seq1 0", e.getMessage().replace('\t', ' '));
    }
  }

  public void testBad3() {
    final String coStr = ""
      + "# a comment" + LS
      + "seq1 0 xxx" + LS
      + "seq1 1 2.0" + LS
      + "seq1 2 0.25" + LS
      + "seq1 5 1.5" + LS
      ;
    try {
      new Corrections(new ByteArrayInputStream(coStr.replace(' ', '\t').getBytes()), 1);
      fail();
    } catch (final Exception e) {
      assertEquals("seq1 0 xxx", e.getMessage().replace('\t', ' '));
    }
  }


}
