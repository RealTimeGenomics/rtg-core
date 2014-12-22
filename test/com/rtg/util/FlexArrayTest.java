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

import junit.framework.TestCase;

/**
 * Test class
 */
public class FlexArrayTest extends TestCase {

  public void testDefault() {
    final FlexArray<Integer> intArr = new FlexArray<>();
    assertEquals(10, intArr.capacity());
  }

  public void testNewSize() {
    assertEquals(2, FlexArray.newSize(0, 1));
    assertEquals(9, FlexArray.newSize(0, 8));
    assertEquals(18, FlexArray.newSize(10, 13));
    assertEquals(Integer.MAX_VALUE, FlexArray.newSize(10, Integer.MAX_VALUE));
  }

  //basic test
  public void test() {
    final FlexArray<Integer> intArr = new FlexArray<>(0);
    assertEquals(0, intArr.capacity());
    intArr.set(5, 7);
    assertEquals(9, intArr.capacity());
    checkSingleSet(intArr, 7, 0);
    checkSingleSet(intArr, 7, 8);
    checkSingleSet(intArr, 7, 15);
  }

  //set 1 so small capacity increase
  public void test1() {
    final FlexArray<Integer> intArr = new FlexArray<>(0);
    assertEquals(0, intArr.capacity());
    intArr.set(5, 1);
    assertEquals(2, intArr.capacity());
    checkSingleSet(intArr, 1, 0);
    checkSingleSet(intArr, 1, 1);
    checkSingleSet(intArr, 1, 15);
  }

  // border case for setting capacity
  public void test2() {
    final FlexArray<Integer> intArr = new FlexArray<>(4);
    assertEquals(4, intArr.capacity());
    intArr.set(5, 4);
    assertEquals(6, intArr.capacity());
    checkSingleSet(intArr, 4, 0);
    checkSingleSet(intArr, 4, 1);
    checkSingleSet(intArr, 4, 15);
  }


  void checkSingleSet(FlexArray<Integer> intArr, int index, int arrLength) {
    final Integer[] ia = intArr.toArray(new Integer[arrLength]);
    assertEquals(index + 1, ia.length);
    assertEquals(5, ia[index].intValue());
    for (int i = 0; i < index + 1; i++) {
      if (i == index) {
        continue;
      }
      assertNull(ia[i]);
    }
    assertEquals(index + 1, intArr.size());
  }
}
