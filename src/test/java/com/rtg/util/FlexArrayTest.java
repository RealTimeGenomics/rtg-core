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
    intArr.set(7, 5);
    assertEquals(9, intArr.capacity());
    checkSingleSet(intArr, 7, 0);
    checkSingleSet(intArr, 7, 8);
    checkSingleSet(intArr, 7, 15);
  }

  //set 1 so small capacity increase
  public void test1() {
    final FlexArray<Integer> intArr = new FlexArray<>(0);
    assertEquals(0, intArr.capacity());
    intArr.set(1, 5);
    assertEquals(2, intArr.capacity());
    checkSingleSet(intArr, 1, 0);
    checkSingleSet(intArr, 1, 1);
    checkSingleSet(intArr, 1, 15);
  }

  // border case for setting capacity
  public void test2() {
    final FlexArray<Integer> intArr = new FlexArray<>(4);
    assertEquals(4, intArr.capacity());
    intArr.set(4, 5);
    assertEquals(6, intArr.capacity());
    checkSingleSet(intArr, 4, 0);
    checkSingleSet(intArr, 4, 1);
    checkSingleSet(intArr, 4, 15);
  }


  void checkSingleSet(FlexArray<Integer> intArr, int index, int arrLength) {
    final Integer[] ia = intArr.toArray(new Integer[arrLength]);
    assertEquals(index + 1, ia.length);
    assertEquals(5, ia[index].intValue());
    for (int i = 0; i < index + 1; ++i) {
      if (i == index) {
        continue;
      }
      assertNull(ia[i]);
    }
    assertEquals(index + 1, intArr.size());
  }
}
