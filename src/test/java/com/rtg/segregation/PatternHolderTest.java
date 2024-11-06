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

package com.rtg.segregation;

import java.util.Arrays;

import com.rtg.util.StringUtils;

import junit.framework.TestCase;

/**
 */
public class PatternHolderTest extends TestCase {

  public void test() {
    //1 10875805  11711404  11011111011 01110000000 X 9 3
    final PatternHolder bedPattern = PatternHolder.fromPatternStrings("11011111011", "01110000000", "X");
    assertEquals(" fa: 11011111011 mo: 01110000000", bedPattern.pattern().toString());
    assertFalse(bedPattern.isNew());
    assertTrue(bedPattern.compatible(bedPattern));
    assertEquals(2, bedPattern.group(0));
    assertEquals(3, bedPattern.group(1));
    assertEquals(1, bedPattern.group(2));
    assertEquals(3, bedPattern.group(3));
    assertEquals(2, bedPattern.group(4));
    assertEquals(2, bedPattern.group(5));
    assertEquals(2, bedPattern.group(6));
    assertEquals(2, bedPattern.group(7));
    assertEquals(0, bedPattern.group(8));
    assertEquals(2, bedPattern.group(9));
    assertEquals(2, bedPattern.group(10));
    final PatternHolder pattern = PatternHolder.fromPatternStrings("11011111011", "01110000100", null);
    assertFalse(bedPattern.compatible(pattern));
    final String[] incompatibleChildren = bedPattern.incompatibleChildren(pattern);
    assertEquals(11, incompatibleChildren.length);
    final String[] expected = StringUtils.split("C,C,C,C,C,C,C,C,I,C,C", ',');
    assertTrue(Arrays.equals(expected, incompatibleChildren));
  }

  public void testNew() {
    final PatternHolder bedPattern = PatternHolder.fromPatternStrings("000???0?0??", "00110001101", "N");
    assertTrue(bedPattern.isNew());
  }
}
