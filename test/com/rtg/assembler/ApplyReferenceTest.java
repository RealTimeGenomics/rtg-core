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
    ApplyReference apply = new ApplyReference(graph, 4, 4, new IntegerOrPercentage(3));
    final List<ApplyReference.AlignmentState> alignments = apply.alignForward(DnaUtils.encodeString("ACCCAGAGACCAGTGTGACCAGTGACGT"), 5, new ContigPosition(1, 5, graph), new GraphTraversions(graph));
    assertEquals(1, alignments.size());
    final ApplyReference.AlignmentState alignmentState = alignments.get(0);
    assertEquals(0, alignmentState.mMismatches.size());
    assertEquals(19, alignmentState.mReferencePosition);
    assertEquals(19, alignmentState.mContigPosition);
    assertEquals(1, alignmentState.mChain.mContigId);

    MemoryPrintStream out = new MemoryPrintStream();
    apply.referenceSequence("foo", DnaUtils.encodeString("ACCCAGAGACCAGTGTGACCAGTGACGT"), out.printStream());
    TestUtils.containsAll(out.toString()
        , "-\tfoo\t19\t0\t0\t19\t0\t19\t1\t[]"
    );
  }

  public void testApplySplit() throws IOException {
    final GraphKmerAttribute graph = GraphMapCliTest.makeGraph(3, new String[]{"ACCCAGAGAC", "GACCAGTGTGACC"}, new long[][]{{1, 2}});
    ApplyReference apply = new ApplyReference(graph, 4, 4, new IntegerOrPercentage(3));
    final List<ApplyReference.AlignmentState> alignments = apply.alignForward(DnaUtils.encodeString("ACCCAGAGACCAGTGTGACCAGTGACGT"), 5, new ContigPosition(1, 5, graph), new GraphTraversions(graph));
    final ApplyReference.AlignmentState alignmentState = alignments.get(0);
    assertEquals(1, alignments.size());
    assertEquals(0, alignmentState.mMismatches.size());
    assertEquals(19, alignmentState.mReferencePosition);
    assertEquals(12, alignmentState.mContigPosition);
    assertEquals(2, alignmentState.mChain.mContigId);
  }
  public void testApplyBranch() throws IOException {
    final GraphKmerAttribute graph = GraphMapCliTest.makeGraph(3, new String[]{"ACCCAGAGAC", "GACCAGTGTGACC", "GACCACCGTAACC"}, new long[][]{{1, 2}, {1, 3}});
    ApplyReference apply = new ApplyReference(graph, 4, 4, new IntegerOrPercentage(3));
    final List<ApplyReference.AlignmentState> alignments = apply.alignForward(DnaUtils.encodeString("ACCCAGAGACCAGTGTGACCAGTGACGT"), 5, new ContigPosition(1, 5, graph), new GraphTraversions(graph));
    final ApplyReference.AlignmentState alignmentState = alignments.get(0);
    assertEquals(1, alignments.size());
    assertEquals(0, alignmentState.mMismatches.size());
    assertEquals(19, alignmentState.mReferencePosition);
    assertEquals(12, alignmentState.mContigPosition);
    assertEquals(2, alignmentState.mChain.mContigId);
  }
  public void testApplyMulti() throws IOException {
    final GraphKmerAttribute graph = GraphMapCliTest.makeGraph(3, new String[]{"ACCCAGAGAC", "GACCAGTGTGACC", "GACCAGGGTCACC"}, new long[][]{{1, 2}, {1, 3}});
    ApplyReference apply = new ApplyReference(graph, 4, 4, new IntegerOrPercentage(3));
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
