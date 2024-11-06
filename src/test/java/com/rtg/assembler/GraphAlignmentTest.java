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

import java.util.Arrays;

import com.rtg.assembler.graph.Graph;

import junit.framework.TestCase;

/**
 */
public class GraphAlignmentTest extends TestCase {
  public void testEquals() {
    final GraphAlignment alignment = new GraphAlignment(Arrays.asList(new AlignmentSection(1, 1, 1), new AlignmentSection(2, 3, 3)), 3, null);
    assertFalse(alignment.equals(null));
    assertFalse(alignment.equals("foo"));
    assertEquals(3, alignment.mScore);
    assertEquals(3, alignment.mScore);
    final GraphAlignment alignment2 = new GraphAlignment(Arrays.asList(new AlignmentSection(1, 1, 4), new AlignmentSection(2, 2, 3)), 3, null);
    assertTrue(alignment.equals(alignment2));
    assertTrue(alignment2.equals(alignment));

    final GraphAlignment alignment3 = new GraphAlignment(Arrays.asList(new AlignmentSection(1, 1, 4), new AlignmentSection(2, 2, 3)), 2, null);
    assertFalse(alignment.equals(alignment3));
    assertFalse(alignment3.equals(alignment));

    final GraphAlignment alignment4 = new GraphAlignment(Arrays.asList(new AlignmentSection(1, 1, 4), new AlignmentSection(2, 2, 4)), 3, null);
    assertFalse(alignment.equals(alignment4));
    assertFalse(alignment4.equals(alignment));

    final GraphAlignment alignment5 = new GraphAlignment(Arrays.asList(new AlignmentSection(1, 2, 4), new AlignmentSection(2, 2, 3)), 3, null);
    assertFalse(alignment.equals(alignment5));
    assertFalse(alignment5.equals(alignment));

    final GraphAlignment alignment6 = new GraphAlignment(Arrays.asList(new AlignmentSection(1, 1, 4), new AlignmentSection(3, 2, 3)), 3, null);
    assertFalse(alignment.equals(alignment6));
    assertFalse(alignment6.equals(alignment));

    final GraphAlignment alignment7 = new GraphAlignment(Arrays.asList(new AlignmentSection(3, 1, 4), new AlignmentSection(2, 2, 3)), 3, null);
    assertFalse(alignment.equals(alignment7));
    assertFalse(alignment7.equals(alignment));

    final GraphAlignment alignment8 = new GraphAlignment(Arrays.asList(new AlignmentSection(1, 1, 4), new AlignmentSection(2, 3, 3), new AlignmentSection(3, 2, 3)), 3, null);
    assertFalse(alignment.equals(alignment8));
    assertFalse(alignment8.equals(alignment));
  }
  public void testToString() {
    final GraphAlignment alignment = new GraphAlignment(Arrays.asList(new AlignmentSection(1, 1, 1), new AlignmentSection(2, 3, 3)), 3, null);
    assertEquals("GraphAlignment{mContigs=[1, 2], mStartPosition=1, mEndPosition=3, mScore=3}", alignment.toString());
  }
  public void testHashCode() {
    final Graph g = GraphMapCliTest.makeGraph(2, new String[] {"ATAT", "CCGG"}, new long[][] {});
    final GraphAlignment alignment = new GraphAlignment(Arrays.asList(new AlignmentSection(1, 1, 1), new AlignmentSection(20000000000L, 3, 3)), 3, g);
    final GraphAlignment alignmentCopy = new GraphAlignment(Arrays.asList(new AlignmentSection(-1, 1, 1), new AlignmentSection(20000000000L, 3, 3)), 3, g);
    assertEquals(alignment.hashCode(), alignmentCopy.hashCode());
    assertEquals(alignment, alignmentCopy);

    final GraphAlignment alignment2 = new GraphAlignment(Arrays.asList(new AlignmentSection(1, 2, 2), new AlignmentSection(20000000000L, 2, 3)), 3, null);
    assertNotSame(alignment.hashCode(), alignment2.hashCode());

    final GraphAlignment alignment3 = new GraphAlignment(Arrays.asList(new AlignmentSection(1, 1, 2), new AlignmentSection(20000000000L, 2, 3)), 3, null);
    assertEquals(alignment.hashCode(), alignment3.hashCode());

    final GraphAlignment alignment4 = new GraphAlignment(Arrays.asList(new AlignmentSection(1, 1, 2), new AlignmentSection(200000000L, 2, 3)), 3, null);
    assertNotSame(alignment.hashCode(), alignment4.hashCode());

    final GraphAlignment alignment5 = new GraphAlignment(Arrays.asList(new AlignmentSection(1, 1, 2), new AlignmentSection(2000000000L, 2, 4)), 3, null);
    assertNotSame(alignment.hashCode(), alignment5.hashCode());
  }

  public void testFromChain() {
    final AlignmentChain forward = new AlignmentChain(2, new AlignmentSection(9, 2, 42), new AlignmentChain(2, new AlignmentSection(2, 2, 2), null));
    final AlignmentChain reverse = new AlignmentChain(4, new AlignmentSection(3, 6, 8), new AlignmentChain(2, new AlignmentSection(2, 2, 2), null));
    final GraphAlignment alignment = new GraphAlignment(forward, reverse, null);
    assertEquals(Arrays.asList(3L, 2L, 9L), alignment.contigs());
    assertEquals(6, alignment.mScore);
    assertEquals(6, alignment.startPosition());
    assertEquals(42, alignment.endPosition());
  }
  public void testFromChainShortForward() {
    final AlignmentChain forward = new AlignmentChain(2, new AlignmentSection(2, 2, 42), null);
    final AlignmentChain reverse = new AlignmentChain(4, new AlignmentSection(3, 6, 8), new AlignmentChain(2, new AlignmentSection(2, 2, 2), null));
    final GraphAlignment alignment = new GraphAlignment(forward, reverse, null);
    assertEquals(Arrays.asList(3L, 2L), alignment.contigs());
    assertEquals(6, alignment.mScore);
    assertEquals(6, alignment.startPosition());
    assertEquals(42, alignment.endPosition());
  }
  public void testFromChainShortReverse() {
    final AlignmentChain forward = new AlignmentChain(2, new AlignmentSection(9, 2, 42), new AlignmentChain(2, new AlignmentSection(2, 2, 2), null));
    final AlignmentChain reverse = new AlignmentChain(4, new AlignmentSection(2, 2, 2), null);
    final GraphAlignment alignment = new GraphAlignment(forward, reverse, null);
    assertEquals(Arrays.asList(2L, 9L), alignment.contigs());
    assertEquals(6, alignment.mScore);
    assertEquals(2, alignment.startPosition());
    assertEquals(42, alignment.endPosition());
  }

  public void testPalindrome() {
    final Graph g = GraphMapCliTest.makeGraph(0, new String[]{"ACGG", "ACGT", "TTTT"}, new long[][]{{1, 2}, {2, 3}});
    final GraphAlignment alignment = new GraphAlignment(Arrays.asList(new AlignmentSection(1, 1, 1), new AlignmentSection(2, 1, 1), new AlignmentSection(3, 1, 1)), 0, g);
    final GraphAlignment alignment2 = new GraphAlignment(Arrays.asList(new AlignmentSection(1, 1, 1), new AlignmentSection(-2, 1, 1), new AlignmentSection(3, 1, 1)), 0, g);
    assertEquals(alignment, alignment2);
    assertEquals(alignment.hashCode(), alignment2.hashCode());
  }
  public void testIsPalindrome() {
    final Graph g = GraphMapCliTest.makeGraph(0, new String[]{"AATT", "ATCT", "AAATT", "AATA", "AAACCCCGGTTT"}, new long[][]{});
    assertTrue(GraphAlignment.isPalindrome(1, g));
    assertFalse(GraphAlignment.isPalindrome(2, g));
    assertFalse(GraphAlignment.isPalindrome(3, g));
    assertFalse(GraphAlignment.isPalindrome(4, g));
    assertFalse(GraphAlignment.isPalindrome(5, g));
  }

  public void testSimpleConstructor() {
    final GraphAlignment g = new GraphAlignment(4, 1, Arrays.asList(1L, 2L, 6L), 40, null);
    assertEquals(4, g.startPosition());
    assertEquals(1, g.endPosition());
    assertEquals(1L, g.startContig());
    assertEquals(6L, g.endContig());
    assertEquals(Arrays.asList(1L, 2L, 6L), g.contigs());

  }
}
