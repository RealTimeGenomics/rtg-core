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


import junit.framework.TestCase;

/**
 */
public class MatchResultTest extends TestCase {

  public final void test() {
    final MatchResult r = new MatchResult(0);
    r.addMatchResult(0, 1, 2, true);
    r.addMatchResult(2, 234, 5, true);
    r.addMatchResult(5, 14, 2, true);
    r.addMatchResult(76, 14, 2, true);
    r.addMatchResult(12, 34, 2, false);
    assertEquals(0, r.getTemplateId(0));
    assertEquals(2, r.getTemplateId(1));
    assertEquals(5, r.getTemplateId(2));
    assertEquals(76, r.getTemplateId(3));
    assertEquals(12, r.getTemplateId(4));
    assertEquals(1, r.getPosition(0));
    assertEquals(34, r.getPosition(4));
    assertEquals(2, r.getEncodedReadId(0));
    assertEquals(2, r.getEncodedReadId(4));
    assertEquals(true, r.isReverse(0));
    assertEquals(false, r.isReverse(4));
    r.sort();
    assertEquals(0, r.getTemplateId(0));
    assertEquals(2, r.getTemplateId(1));
    assertEquals(5, r.getTemplateId(2));
    assertEquals(12, r.getTemplateId(3));
    assertEquals(76, r.getTemplateId(4));
    assertEquals(1, r.getPosition(0));
    assertEquals(14, r.getPosition(4));
    assertEquals(2, r.getEncodedReadId(0));
    assertEquals(2, r.getEncodedReadId(4));
    assertEquals(true, r.isReverse(0));
    assertEquals(true, r.isReverse(4));

  }

}
