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
package com.rtg.variant.match;

import com.rtg.reader.FastaUtils;

import junit.framework.TestCase;

/**
 */
public class MatchTest extends TestCase {

  public void testToString1() {
    final Match ma = new AlignmentMatch(null, null, "ACGT", null, 10, 0, 4, 10, false, false);
    assertEquals("~ACGT~", ma.toString());
  }

  public void testToString2() {
    final Match ma = new AlignmentMatch(null, null, "ACGT", null, 10, 0, 4, 10, true, false);
    assertEquals("ACGT~", ma.toString());
  }

  public void testToString3() {
    final Match ma = new AlignmentMatch(null, null, "ACGT", null, 10, 0, 4, 10, false, true);
    assertEquals("~ACGT", ma.toString());
  }

  public void testToString4() {
    final Match ma = new AlignmentMatch(null, null, "ACGT", null, 10, 0, 4, 10, true, true);
    assertEquals("ACGT", ma.toString());
  }

  public void testQualityString() {
    final Match ma = new AlignmentMatch(null, null, "ACGT", FastaUtils.asciiToRawQuality("!%$#"), 10, 0, 4, 10, false, false);
    assertEquals("[0.00, 0.40, 0.30, 0.20, ]", ma.qualityString());
  }

  public void testCorrection() {
    // see diploidcomplexscorertest.xls for details
    final Match m1 = new AlignmentMatch(null, "ACGT", null, 20, 0, 4, 30);
    assertEquals(0.01099, m1.correction(), 0.00001);
    final Match m2 = new AlignmentMatch(null, "ACGT", null, 20, 0, 4, 10);
    assertEquals(0.109, m2.correction(), 0.00001);
    final Match m3 = new AlignmentMatch(null, "ACGT", "$$#!", 0, 0, 4, 10);
    assertEquals(0.69250, m3.correction(), 0.00001);
    final Match m4 = new AlignmentMatch(null, "", "", 0, 0, 0, 10);
    assertEquals(0.1, m4.correction(), 0.00001);
  }
}

