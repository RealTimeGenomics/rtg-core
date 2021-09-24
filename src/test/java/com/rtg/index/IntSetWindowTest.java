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

import static com.rtg.util.StringUtils.LS;

import java.io.IOException;

import junit.framework.TestCase;



/**
 */
public class IntSetWindowTest extends TestCase {

  public void test() throws IOException {
    final StringBuilder sb = new StringBuilder();
    final IntSetCaller caller = new IntSetCaller() {
      @Override
      public void call(final int v) {
        sb.append(v).append(" ");
      }
    };
    final IntSet is = new IntSetWindow(10, 3, 1, caller);
    is.globalIntegrity();
    assertEquals("IntSetWindow 1" + LS + "0: " + LS, is.toString());
    is.add(9);
    is.add(0);
    is.add(1);
    is.add(9);
    is.add(1);
    is.globalIntegrity();
    assertEquals("IntSetWindow 1" + LS + "3: 9 0 1 " + LS, is.toString());
    is.iterateClear();
    is.globalIntegrity();
    assertEquals("IntSetWindow 1" + LS + "0: " + LS, is.toString());
    assertEquals("9 0 1 ", sb.toString());

    is.add(0);
    is.add(1);
    is.add(1);
    is.globalIntegrity();
    assertEquals("IntSetWindow 1" + LS + "2: 0 1 " + LS, is.toString());
  }

  public void test3() throws IOException {
    final StringBuilder sb = new StringBuilder();
    final IntSetCaller caller = new IntSetCaller() {
      @Override
      public void call(final int v) {
        sb.append(v).append(" ");
      }
    };
    final IntSet is = new IntSetWindow(10, 3, 3, caller);
    is.globalIntegrity();
    final String exp0 = ""
      + "IntSetWindow 3" + LS
      + "0: " + LS
      + "0: " + LS
      + "0: " + LS
      ;
    assertEquals(exp0, is.toString());
    is.add(9);
    is.add(0);
    is.add(1);
    is.add(9);
    is.add(1);
    is.globalIntegrity();
    final String exp1a = ""
      + "IntSetWindow 3" + LS
      + "0: " + LS
      + "0: " + LS
      + "3: 9 0 1 " + LS
      ;
    assertEquals(exp1a, is.toString());

    is.iterateClear();
    is.globalIntegrity();
    final String exp1b = ""
      + "IntSetWindow 3" + LS
      + "0: " + LS
      + "3: 9 0 1 " + LS
      + "0: " + LS
      ;
    assertEquals(exp1b, is.toString());
    assertEquals("", sb.toString());

    is.iterateClear();
    is.add(5);
    is.globalIntegrity();
    final String exp1c = ""
      + "IntSetWindow 3" + LS
      + "3: 9 0 1 " + LS
      + "0: " + LS
      + "1: 5 " + LS
      ;
    assertEquals(exp1c, is.toString());
    assertEquals("", sb.toString());

    is.iterateClear();
    is.globalIntegrity();
    final String exp1d = ""
      + "IntSetWindow 3" + LS
      + "0: " + LS
      + "1: 5 " + LS
      + "0: " + LS
      ;
    assertEquals(exp1d, is.toString());
    assertEquals("9 0 1 ", sb.toString());

    is.add(0);
    is.add(1);
    is.add(1);
    is.globalIntegrity();
    final String exp2 = ""
      + "IntSetWindow 3" + LS
      + "0: " + LS
      + "1: 5 " + LS
      + "2: 0 1 " + LS
      ;
    assertEquals(exp2, is.toString());

    is.iterateClearAll();
    is.globalIntegrity();
    final String exp2b = ""
      + "IntSetWindow 3" + LS
      + "0: " + LS
      + "0: " + LS
      + "0: " + LS
      ;
    assertEquals(exp2b, is.toString());
    assertEquals("9 0 1 5 0 1 ", sb.toString());
  }

  public void testEmpty() throws IOException {
    final StringBuilder sb = new StringBuilder();
    final IntSetCaller caller = new IntSetCaller() {
      @Override
      public void call(final int v) {
        sb.append(v).append(" ");
      }
    };
    final IntSet is = new IntSetWindow(10, 0, 1, caller);
    is.globalIntegrity();
    assertEquals("IntSetWindow 1" + LS + "0: " + LS, is.toString());
    is.add(9);
    is.add(0);
    is.add(1);
    is.add(9);
    is.add(1);
    is.globalIntegrity();
    assertEquals("IntSetWindow 1" + LS + "0: " + LS, is.toString());
    is.iterateClear();
    is.globalIntegrity();
    assertEquals("IntSetWindow 1" + LS + "0: " + LS, is.toString());
    assertEquals("9 0 1 9 1 ", sb.toString());
  }

  public void testBad() throws IOException {
    final IntSetCaller caller = new IntSetCaller() {
      @Override
      public void call(final int v) {
        //do nothing
      }
    };
    final IntSet is = new IntSetWindow(10, 3, 1, caller);
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
      new IntSetWindow(Integer.MAX_VALUE + 1L, 3, 1, caller);
      fail();
    } catch (final RuntimeException e) {
      assertEquals("Range too large:2147483648", e.getMessage());
    }

  }
}
