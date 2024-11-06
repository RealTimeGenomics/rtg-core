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
package com.rtg.ngs;

import com.rtg.util.TestUtils;

import junit.framework.TestCase;

/**
 */
public class NgsMaskParamsGeneralTest extends TestCase {

  NgsMaskParamsGeneral getParams(final int wordsize, final int substitutions, final int indels, final int indelLength) {
    return new NgsMaskParamsGeneral(wordsize, substitutions, indels, indelLength);
  }

  public void testEquals() {
    final int readLength = 36;
    final NgsMaskParams a1 = getParams(12, 2, 1, 1);
    assertTrue(a1.isValid(readLength));
    final NgsMaskParams a2 = getParams(12, 2, 1, 1);
    assertTrue(a2.isValid(readLength));
    final NgsMaskParams z = getParams(12, 2, 2, 2);
    assertTrue(z.isValid(readLength));
    final NgsMaskParams b = getParams(12, 2, 2, 1);
    assertTrue(b.isValid(readLength));
    final NgsMaskParams c = getParams(12, 3, 1, 1);
    assertTrue(c.isValid(readLength));
    final NgsMaskParams d = getParams(13, 2, 1, 1);
    assertTrue(d.isValid(readLength));
    TestUtils.equalsHashTest(new NgsMaskParams[][] {{a1, a2}, {z}, {b}, {c}, {d}});
  }


  public void test() throws Exception {
    final NgsMaskParamsGeneral a = getParams(12, 3, 2, 1);
    a.integrity();
    assertEquals("General Mask: w=12 s=3 i=2 l=1", a.toString());
    assertNotNull(a.maskFactory(36));
    a.close();
  }

  public void testInvalid() {
    final NgsMaskParams a = getParams(32, 1, 1, 1);
    assertFalse(a.isValid(33));
  }
}
