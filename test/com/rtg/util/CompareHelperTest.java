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

package com.rtg.util;

import java.util.Arrays;
import java.util.List;

import junit.framework.TestCase;

/**
 */
public class CompareHelperTest extends TestCase {
  public void testDefault() {
    assertEquals(0, new CompareHelper().result());
  }

  public void testChaining() {
    CompareHelper ch = new CompareHelper();
    assertEquals(ch, ch.compare(3, 3));
    assertEquals(0, ch.result());
    assertEquals(ch, ch.compare("Foo", "Foo")
        .compare(5.0, 8.0)
        .compare(5.0, 5.0));
    assertEquals(Double.compare(5.0, 8.0),  ch.result());
  }

  public void testString() {
    assertEquals("foo".compareTo("Foo")  , new CompareHelper()
        .compare(30, 30)
        .compare("foo", "Foo")
        .compare(10.0, 50.0)
        .result()
    );
  }
  public void testListCompare() {
    List<Integer> first = Arrays.asList(1, 2, 3, 4);
    List<Integer> second = Arrays.asList(1, 2, 3, 4, 5);
    List<Integer> third = Arrays.asList(1, 2, 4, 4, 5);
    List<Integer> fourth = Arrays.asList(-1, 2, 4, 4, 5);
    assertEquals(0, new CompareHelper().compareList(first, first).result());
    assertEquals(0, new CompareHelper().compareList(second, second).result());
    assertEquals(Integer.valueOf(4).compareTo(5), new CompareHelper().compareList(first, second).result());
    assertEquals(Integer.valueOf(3).compareTo(4), new CompareHelper().compareList(first, third).result());
    assertEquals(Integer.valueOf(1).compareTo(-1), new CompareHelper().compareList(first, fourth).result());
  }
  public void testListCompareAlreadyDifferent() {
    List<Integer> first = Arrays.asList(3, 2, 3, 4);
    List<Integer> second = Arrays.asList(1, 2, 3, 4);
    assertEquals("foo".compareTo("bar"), new CompareHelper().compare("foo", "bar").compareList(first, second).result());
  }
}
