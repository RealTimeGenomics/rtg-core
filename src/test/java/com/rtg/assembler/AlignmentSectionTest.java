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
