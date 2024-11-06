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
 */
public class MemoTest extends TestCase {

  public void test() {
    final Function<Integer> function = new Function<Integer>() {
      @Override
      public Integer fn(int i) {
        return i + 42;
      }
    };
    final Memo<Integer> im = new Memo<>(function);
    im.globalIntegrity();
    assertEquals("Memo empty" + LS, im.toString());
    assertEquals(Integer.valueOf(42), im.fn(0));
    im.globalIntegrity();
    assertEquals(Integer.valueOf(44), im.fn(2));
    im.globalIntegrity();
    assertEquals(Integer.valueOf(37), im.fn(-5));
    im.globalIntegrity();
    assertEquals("Memo [-5..2]" + LS + "37, 38, 39, 40, 41, 42, 43, 44" + LS, im.toString());

    assertEquals(Integer.valueOf(35), im.fn(-7));
    im.globalIntegrity();
    assertEquals(Integer.valueOf(46), im.fn(4));
    im.globalIntegrity();
    assertEquals("Memo [-7..4]" + LS + "35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46" + LS, im.toString());
    assertEquals(Integer.valueOf(42), im.fn(0));
    assertEquals(Integer.valueOf(44), im.fn(2));
    assertEquals(Integer.valueOf(37), im.fn(-5));
  }
}
