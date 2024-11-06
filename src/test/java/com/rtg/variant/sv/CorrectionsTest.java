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
