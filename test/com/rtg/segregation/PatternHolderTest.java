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

import java.util.Arrays;

import com.rtg.util.StringUtils;

import junit.framework.TestCase;

/**
 */
public class PatternHolderTest extends TestCase {

  public void test() throws MismatchingPloidyException {
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

  public void testNew() throws MismatchingPloidyException {
    final PatternHolder bedPattern = PatternHolder.fromPatternStrings("000???0?0??", "00110001101", "N");
    assertTrue(bedPattern.isNew());
  }
}
