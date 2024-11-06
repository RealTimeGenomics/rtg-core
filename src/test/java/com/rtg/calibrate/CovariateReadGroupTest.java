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

package com.rtg.calibrate;

import htsjdk.samtools.SAMRecord;

import junit.framework.TestCase;

/**
 */
public class CovariateReadGroupTest extends TestCase {

  public void test() {
    final CovariateReadGroup rg = new CovariateReadGroup();
    assertEquals("readgroup", rg.name());
    assertEquals(CovariateEnum.READGROUP, rg.getType());
    assertEquals(0, rg.size());
    assertEquals(0, rg.parse("foo"));
    assertEquals(1, rg.newSize());
    assertEquals(0, rg.size());
    assertEquals(true, rg.sizeChanged());
    rg.resized();
    assertEquals(false, rg.sizeChanged());
    assertEquals(1, rg.size());
    assertEquals(1, rg.parse("bar"));
    assertEquals(2, rg.newSize());
    assertEquals(0, rg.parse("foo"));
    assertEquals(2, rg.newSize());
    assertEquals("foo", rg.valueString(0));
    assertEquals("bar", rg.valueString(1));
    assertEquals(-1, rg.valueOf("dasfdfas"));
  }


  public void testValue() {
    final CovariateReadGroup rg = new CovariateReadGroup();
    final SAMRecord sam = new SAMRecord(null);
    assertEquals(0, rg.value(sam, null));
    assertEquals(1, rg.newSize());
    assertEquals("unspecified", rg.valueString(0));

    sam.setAttribute("RG", "foo");
    assertEquals(1, rg.value(sam, null));
    assertEquals(2, rg.newSize());
    assertEquals("foo", rg.valueString(1));
  }
}
