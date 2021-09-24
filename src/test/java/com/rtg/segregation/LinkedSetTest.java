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

package com.rtg.segregation;

import java.util.Iterator;

import junit.framework.TestCase;

/**
 */
public class LinkedSetTest extends TestCase {

  public void test() {
    final LinkedSet<String> ls = new LinkedSet<>();
    ls.globalIntegrity();
    assertTrue(ls.add("foo"));
    ls.globalIntegrity();
    assertTrue(ls.add("baz"));
    ls.globalIntegrity();
    assertTrue(ls.add("bar"));
    ls.globalIntegrity();
    assertFalse(ls.add("foo"));
    ls.globalIntegrity();

    check(ls.iterator("foo", "bar"), "baz bar ");
    check(ls.iterator("foo", "foo"), "");
    check(ls.iterator("baz", "bar"), "bar ");
  }

  private void check(Iterator<String> it, final String exp) {
    final StringBuilder sb = new StringBuilder();
    while (it.hasNext()) {
      sb.append(it.next()).append(" ");
    }
    assertEquals(exp, sb.toString());
  }
}
