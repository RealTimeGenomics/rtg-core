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

package com.rtg.assembler;

import junit.framework.TestCase;

/**
 */
public class AlignmentSectionTest extends TestCase {
  public void testToString() {
    assertEquals("AlignmentSection{mContig=1, mStartPosition=2, mEndPosition=3}", new AlignmentSection(1, 2, 3).toString());
  }

  public void testEquals() {
    final AlignmentSection alignmentSection = new AlignmentSection(1, 2, 3);
    assertFalse(alignmentSection.equals(null));
    assertFalse(alignmentSection.equals("foo"));
    assertTrue(alignmentSection.equals(alignmentSection));

    final AlignmentSection alignmentSection2 = new AlignmentSection(2, 2, 3);
    assertFalse(alignmentSection.equals(alignmentSection2));
    assertFalse(alignmentSection2.equals(alignmentSection));
    final AlignmentSection alignmentSection3 = new AlignmentSection(1, 3, 3);
    assertFalse(alignmentSection.equals(alignmentSection3));
    assertFalse(alignmentSection3.equals(alignmentSection));
    final AlignmentSection alignmentSection4 = new AlignmentSection(1, 2, 2);
    assertFalse(alignmentSection.equals(alignmentSection4));
    assertFalse(alignmentSection4.equals(alignmentSection));
    final AlignmentSection alignmentSection5 = new AlignmentSection(1, 2, 3);
    assertTrue(alignmentSection.equals(alignmentSection5));
  }

  public void testHashCode() {
    final AlignmentSection alignmentSection = new AlignmentSection(20000000000L, 2, 3);
    final AlignmentSection alignmentSection2 = new AlignmentSection(2, 2, 3);
    assertTrue(alignmentSection.hashCode() != alignmentSection2.hashCode());
    assertEquals(21354309, alignmentSection.hashCode());

  }
}
