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
package com.rtg.util.memo;

import static com.rtg.util.StringUtils.LS;

import junit.framework.TestCase;

/**
 *
 *
 */
public class DoubleMemoTest extends TestCase {

  public void test() {
    final DoubleFunction function = new DoubleFunction() {
      @Override
      public double fn(int i) {
        return i + 42.3;
      }
    };
    final DoubleMemo im = new DoubleMemo(function);
    im.globalIntegrity();
    assertEquals("IntMemo empty" + LS, im.toString());
    assertEquals(42.3, im.fn(0));
    im.globalIntegrity();
    assertEquals(44.3, im.fn(2));
    im.globalIntegrity();
    assertEquals(37.3, im.fn(-5));
    im.globalIntegrity();
    assertEquals("IntMemo [-5..2]" + LS + "37.3, 38.3, 39.3, 40.3, 41.3, 42.3, 43.3, 44.3" + LS, im.toString());

    assertEquals(35.3, im.fn(-7));
    im.globalIntegrity();
    assertEquals(46.3, im.fn(4));
    im.globalIntegrity();
    assertEquals("IntMemo [-7..4]" + LS + "35.3, 36.3, 37.3, 38.3, 39.3, 40.3, 41.3, 42.3, 43.3, 44.3, 45.3, 46.3" + LS, im.toString());
    assertEquals(42.3, im.fn(0));
    assertEquals(44.3, im.fn(2));
    assertEquals(37.3, im.fn(-5));
  }
}
