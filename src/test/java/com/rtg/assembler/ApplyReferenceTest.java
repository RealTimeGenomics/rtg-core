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

import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import com.rtg.assembler.graph.implementation.GraphKmerAttribute;
import com.rtg.mode.DnaUtils;
import com.rtg.util.IntegerOrPercentage;
import com.rtg.util.TestUtils;
import com.rtg.util.io.MemoryPrintStream;

import junit.framework.TestCase;

/**
 */
public class ApplyReferenceTest extends TestCase {
  public void testApply() throws IOException {
    final GraphKmerAttribute graph = GraphMapCliTest.makeGraph(3, new String[]{"ACCCAGAGACCAGTGTGACC"}, new long[][]{});
    final ApplyReference apply = new ApplyReference(graph, 4, 4, new IntegerOrPercentage(3));
    final List<ApplyReference.AlignmentState> alignments = apply.alignForward(DnaUtils.encodeString("ACCCAGAGACCAGTGTGACCAGTGACGT"), 5, new ContigPosition(1, 5, graph), new GraphTraversions(graph));
    assertEquals(1, alignments.size());
    final ApplyReference.AlignmentState alignmentState = alignments.get(0);
    assertEquals(0, alignmentState.mMismatches.size());
    assertEquals(19, alignmentState.mReferencePosition);
    assertEquals(19, alignmentState.mContigPosition);
    assertEquals(1, alignmentState.mChain.mContigId);

    final MemoryPrintStream out = new MemoryPrintStream();
    apply.referenceSequence("foo", DnaUtils.encodeString("ACCCAGAGACCAGTGTGACCAGTGACGT"), out.printStream());
    TestUtils.containsAll(out.toString()
        , "-\tfoo\t19\t0\t0\t19\t0\t19\t1\t[]"
    );
  }

  public void testApplySplit() {
    final GraphKmerAttribute graph = GraphMapCliTest.makeGraph(3, new String[]{"ACCCAGAGAC", "GACCAGTGTGACC"}, new long[][]{{1, 2}});
    final ApplyReference apply = new ApplyReference(graph, 4, 4, new IntegerOrPercentage(3));
    final List<ApplyReference.AlignmentState> alignments = apply.alignForward(DnaUtils.encodeString("ACCCAGAGACCAGTGTGACCAGTGACGT"), 5, new ContigPosition(1, 5, graph), new GraphTraversions(graph));
    final ApplyReference.AlignmentState alignmentState = alignments.get(0);
    assertEquals(1, alignments.size());
    assertEquals(0, alignmentState.mMismatches.size());
    assertEquals(19, alignmentState.mReferencePosition);
    assertEquals(12, alignmentState.mContigPosition);
    assertEquals(2, alignmentState.mChain.mContigId);
  }
  public void testApplyBranch() {
    final GraphKmerAttribute graph = GraphMapCliTest.makeGraph(3, new String[]{"ACCCAGAGAC", "GACCAGTGTGACC", "GACCACCGTAACC"}, new long[][]{{1, 2}, {1, 3}});
    final ApplyReference apply = new ApplyReference(graph, 4, 4, new IntegerOrPercentage(3));
    final List<ApplyReference.AlignmentState> alignments = apply.alignForward(DnaUtils.encodeString("ACCCAGAGACCAGTGTGACCAGTGACGT"), 5, new ContigPosition(1, 5, graph), new GraphTraversions(graph));
    final ApplyReference.AlignmentState alignmentState = alignments.get(0);
    assertEquals(1, alignments.size());
    assertEquals(0, alignmentState.mMismatches.size());
    assertEquals(19, alignmentState.mReferencePosition);
    assertEquals(12, alignmentState.mContigPosition);
    assertEquals(2, alignmentState.mChain.mContigId);
  }
  public void testApplyMulti() {
    final GraphKmerAttribute graph = GraphMapCliTest.makeGraph(3, new String[]{"ACCCAGAGAC", "GACCAGTGTGACC", "GACCAGGGTCACC"}, new long[][]{{1, 2}, {1, 3}});
    final ApplyReference apply = new ApplyReference(graph, 4, 4, new IntegerOrPercentage(3));
    final List<ApplyReference.AlignmentState> alignments = apply.alignForward(DnaUtils.encodeString("ACCCAGAGACCAGTGTGACCAGTGACGT"), 5, new ContigPosition(1, 5, graph), new GraphTraversions(graph));
    assertEquals(2, alignments.size());
    ApplyReference.AlignmentState alignmentState = alignments.get(0);
    assertEquals(Collections.<Integer>emptyList(), alignmentState.mMismatches);
    assertEquals(19, alignmentState.mReferencePosition);
    assertEquals(12, alignmentState.mContigPosition);
    assertEquals(2, alignmentState.mChain.mContigId);

    alignmentState = alignments.get(1);
    assertEquals(Arrays.asList(13, 16), alignmentState.mMismatches);
    assertEquals(19, alignmentState.mReferencePosition);
    assertEquals(12, alignmentState.mContigPosition);
    assertEquals(3, alignmentState.mChain.mContigId);
  }
}
