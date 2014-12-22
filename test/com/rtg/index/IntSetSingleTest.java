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
