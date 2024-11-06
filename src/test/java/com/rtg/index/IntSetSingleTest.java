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
package com.rtg.index;


import java.io.IOException;

import com.rtg.util.StringUtils;

import junit.framework.TestCase;



/**
 */
public class IntSetSingleTest extends TestCase {

  public void test1() throws IOException {
    final StringBuilder sb = new StringBuilder();
    final IntSetCaller caller = new IntSetCaller() {
      @Override
      public void call(final int v) {
        sb.append(v).append(" ");
      }
    };
    final IntSetSingle is = new IntSetSingle(10, 3, caller);
    is.globalIntegrity();
    assertEquals("IntSet 0" + StringUtils.LS + StringUtils.LS, is.toString());
    is.add(9);
    is.add(0);
    is.add(1);
    is.add(9);
    is.add(1);
    is.globalIntegrity();
    assertEquals("IntSet 3" + StringUtils.LS + "9 0 1 " + StringUtils.LS, is.toString());
    is.iterateClear();
    is.globalIntegrity();
    assertEquals("IntSet 0" + StringUtils.LS + StringUtils.LS, is.toString());
    assertEquals("9 0 1 ", sb.toString());

    is.add(0);
    is.add(1);
    is.add(1);
    is.globalIntegrity();
    assertEquals("IntSet 2" + StringUtils.LS + "0 1 " + StringUtils.LS, is.toString());
  }

  public void test2() throws IOException {
    final StringBuilder sb = new StringBuilder();
    final IntSetCaller caller = new IntSetCaller() {
      @Override
      public void call(final int v) {
        sb.append(v).append(" ");
      }
    };
    final IntSetSingle is = new IntSetSingle(10, 3, caller);
    is.globalIntegrity();
    assertEquals("IntSet 0" + StringUtils.LS + StringUtils.LS, is.toString());
    is.add(9);
    is.add(0);
    is.add(1);
    is.add(9);
    is.add(1);
    is.globalIntegrity();
    assertEquals("IntSet 3" + StringUtils.LS + "9 0 1 " + StringUtils.LS, is.toString());
    is.iterateClearAll();
    is.globalIntegrity();
    assertEquals("IntSet 0" + StringUtils.LS + StringUtils.LS, is.toString());
    assertEquals("9 0 1 ", sb.toString());

    is.add(0);
    is.add(1);
    is.add(1);
    is.globalIntegrity();
    assertEquals("IntSet 2" + StringUtils.LS + "0 1 " + StringUtils.LS, is.toString());
  }

  public void testEmpty() throws IOException {
    final StringBuilder sb = new StringBuilder();
    final IntSetCaller caller = new IntSetCaller() {
      @Override
      public void call(final int v) {
        sb.append(v).append(" ");
      }
    };
    final IntSetSingle is = new IntSetSingle(10, 0, caller);
    is.globalIntegrity();
    assertEquals("IntSet 0" + StringUtils.LS + StringUtils.LS, is.toString());
    is.add(9);
    is.add(0);
    is.add(1);
    is.add(9);
    is.add(1);
    is.globalIntegrity();
    assertEquals("IntSet 0" + StringUtils.LS + StringUtils.LS, is.toString());
    is.iterateClear();
    is.globalIntegrity();
    assertEquals("IntSet 0" + StringUtils.LS + StringUtils.LS, is.toString());
    assertEquals("9 0 1 9 1 ", sb.toString());
  }

  public void testBad() throws IOException {
    final IntSetCaller caller = new IntSetCaller() {
      @Override
      public void call(final int v) {
        //do nothing
      }
    };
    final IntSetSingle is = new IntSetSingle(10, 3, caller);
    try {
      is.add(-1);
      fail();
    } catch (final RuntimeException e) {
      //
    }
    try {
      is.add(10);
      fail();
    } catch (final RuntimeException e) {
      //
    }

    try {
      new IntSetSingle(Integer.MAX_VALUE + 1L, 3, caller);
      fail();
    } catch (final RuntimeException e) {
      assertEquals("Range too large:2147483648", e.getMessage());
    }
  }
}
