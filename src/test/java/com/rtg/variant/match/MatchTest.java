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
    final Match m4 = new AlignmentMatch(null, "", "", 10, 0, 0, 10);
    assertEquals(0.1, m4.correction(), 0.19);
  }
}

