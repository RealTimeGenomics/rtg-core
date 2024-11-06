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
package com.rtg.simulation;

import junit.framework.TestCase;

/**
 */
public class ReadSimEvalStatisticsTest extends TestCase {


  public void test() {
    final ReadSimEvalStatistics s = new ReadSimEvalStatistics(20);
    assertEquals(s.length(), 20);
    s.found(0);
    assertTrue(s.isFound(0));
    s.mated(0);

    assertTrue(s.isMated(0));

    assertFalse(s.isFound(1));
    //assertFalse(s.isBetter(1));
    assertFalse(s.isMultiple(1));
    assertFalse(s.isUnmapped(1));
    assertFalse(s.isUnmated(1));

    s.multiple(1);
    assertTrue(s.isMultiple(1));
    s.unmapped(2);
    assertTrue(s.isUnmapped(2));
    s.unmated(3);
    assertTrue(s.isUnmated(3));

    assertFalse(s.isMapped(4));
    s.mapped(4);
    assertTrue(s.isMapped(4));

  }

}
