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

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;

import com.rtg.assembler.graph.Graph;
import com.rtg.assembler.graph.io.GraphReader;
import com.rtg.index.hash.ExactHashFunction;
import com.rtg.mode.DNA;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.util.IntegerOrPercentage;
import com.rtg.util.StringUtils;
import com.rtg.util.store.StoreDirProxy;

/**
 */
public class ApplyReference {
  final GraphIndex mIndex;
  final Graph mGraph;
  final GraphTraversions mTraversions;
  final IntegerOrPercentage mMismatches;
  ApplyReference(Graph graph, int word, int step, IntegerOrPercentage mismatches) {
    mIndex = new GraphIndex(graph, word, step);
    mGraph = graph;
    mTraversions = new GraphTraversions(mGraph);
    mMismatches = mismatches;
  }

  void referenceSequence(String referenceName, byte[] reference, PrintStream out) throws IOException {
    int referencePosition = 0;
    final ExactHashFunction search = mIndex.getSearchFunction();
    while (referencePosition < reference.length) {
      final byte base = reference[referencePosition];
      final byte code = (byte) (base - 1);
      if (code < 0) {
        search.reset();
      } else {
        search.hashStep(code);
      }
      int nextRefPos = referencePosition;
      String initialPrefix = "-";
      for (ContigPosition contigPosition : mIndex.hashHits(search, mGraph)) {
        final List<AlignmentState> alignmentStates = alignForward(reference, referencePosition, contigPosition, mTraversions);
        if (alignmentStates.size() > 0) {
          out.print(output(initialPrefix, search.getWindowSize(), referenceName, referencePosition, contigPosition, alignmentStates));
          nextRefPos = Math.max(alignmentStates.get(0).mReferencePosition, nextRefPos);
          search.reset();
        }
        initialPrefix = "\\";
      }
      referencePosition = nextRefPos;
      referencePosition++;
    }
  }
  void tabbed(StringBuilder sb, String... vals) {
    String join = "";
    for (String val : vals) {
      sb.append(join).append(val);
      join = "\t";
    }
  }

  String output(String initialPrefix, int windowSize, String referenceName, int refPos, ContigPosition contigPos, List<AlignmentState> states) {
    final StringBuilder sb = new StringBuilder();
    String prefix = initialPrefix;
    for (AlignmentState state : states) {
      final int start = refPos - windowSize + 1;
      final int end = states.get(0).mReferencePosition;
      final int length = states.get(0).mReferencePosition - start;

      tabbed(sb, prefix, referenceName
          , "" + length
          , "" + state.mMismatchCount
          , "" + start
          , "" + end
          , "" + (contigPos.mPosition - windowSize + 1)
          , "" + state.mContigPosition
          , "" + state.mChain
          , "" + state.mMismatches);
      prefix = "\\";
      sb.append(StringUtils.LS);
    }
    return sb.toString();

  }
  static class ReferenceChain {
    long mContigId;
    ReferenceChain mPrevious;
    ReferenceChain(long contigId, ReferenceChain previous) {
      mContigId = contigId;
      mPrevious = previous;
    }
    public String toString() {
      final  StringBuilder sb = new StringBuilder();
      ReferenceChain current = this;
      String join = "";
      while (current != null) {
        sb.insert(0, join);
        sb.insert(0, current.mContigId);
        join = ",";
        current = current.mPrevious;
      }
      return sb.toString();
    }

  }
  static class MismatchComparator implements Comparator<AlignmentState>, Serializable {
    @Override
    public int compare(AlignmentState o1, AlignmentState o2) {
      return Integer.valueOf(o1.numberMismatches()).compareTo(o2.numberMismatches());
    }
  }
  static class AlignmentState {
    int mReferencePosition;
    int mReferenceOrigin;
    int mContigPosition;
    ReferenceChain mChain;
    int mMismatchCount;
    Queue<Integer> mMismatches = new LinkedList<>();

    long contigId() {
      return mChain.mContigId;
    }

    public AlignmentState(int refPos, ReferenceChain chain, int position, int referenceOrigin, int mismatchCount) {
      mChain = chain;
      mReferencePosition = refPos;
      mContigPosition = position;
      mReferenceOrigin = referenceOrigin;
      mMismatchCount = mismatchCount;
    }
    int numberMismatches() {
      return mMismatches.size() + mMismatchCount;
    }
    void addMismatches(Collection<Integer> mismatches) {
      mMismatches.addAll(mismatches);
      trimMismatches();
    }

    void trimMismatches() {
      Integer peeked;
      while ((peeked = mMismatches.peek()) != null) {
        if (mReferencePosition - peeked > 100) {
          mMismatches.remove();
          mMismatchCount++;
        } else {
          break;
        }
      }
    }

    boolean isValid(int mismatches) {
      return mMismatches.size() < mismatches;
    }

    @Override
    public String toString() {
      return "AlignmentState{"
          + "mReferenceOrigin=" + mReferenceOrigin
          + ", mReferencePosition=" + mReferencePosition
          + ", mContigPosition=" + mContigPosition
          + ", mChain=" + mChain
          + ", mMismatchCount=" + mMismatchCount
          + ", mMismatches=" + mMismatches + '}';
    }

    List<AlignmentState> step(byte[] reference, Graph graph, GraphTraversions traversions) {
      final int nextContigPos = mContigPosition + 1;
      final int nextRefPos = mReferencePosition + 1;
      if (nextRefPos >= reference.length) {
        return Collections.emptyList();
      }

      final long contigId = contigId();
      final List<Integer> nextMismatches;
      if (reference[mReferencePosition] != DNA.N.ordinal() && reference[mReferencePosition] != graph.nt(contigId, mContigPosition)) {
        nextMismatches = Collections.singletonList(mReferencePosition);
      } else {
        nextMismatches = Collections.emptyList();
      }

      if (nextContigPos < graph.contigLength(contigId)) {
        final AlignmentState nextState = new AlignmentState(nextRefPos, mChain, nextContigPos, mReferenceOrigin, mMismatchCount);
        nextState.addMismatches(mMismatches);
        nextState.addMismatches(nextMismatches);
        return Collections.singletonList(nextState);
      } else {
        final List<AlignmentState> nextStates = new ArrayList<>();
        for (Long nextContig : traversions.get(contigId).next()) {
          final AlignmentState state = new AlignmentState(nextRefPos, new ReferenceChain(nextContig, mChain), graph.contigOverlap(), mReferenceOrigin, mMismatchCount);
          state.addMismatches(mMismatches);
          state.addMismatches(nextMismatches);
          nextStates.add(state);
        }
        return nextStates;
      }
    }
  }

  List<AlignmentState> alignForward(byte[] reference, int refStart, ContigPosition startPos, GraphTraversions traversion) {
    List<AlignmentState> states = new ArrayList<>();
    ArrayList<AlignmentState> finalStates = new ArrayList<>();
    states.add(new AlignmentState(refStart, new ReferenceChain(startPos.mContigId, null), startPos.mPosition, refStart, 0));
    while (states.size() > 0) {
      finalStates = new ArrayList<>();
      final ArrayList<AlignmentState> nextStates = new ArrayList<>();
      for (AlignmentState state : states) {
        final List<AlignmentState> steppedStates = state.step(reference, mGraph, traversion);
        if (steppedStates.size() == 0) {
          finalStates.add(state);
        } else {
          for (AlignmentState step : steppedStates) {
            step.trimMismatches();
            if (step.isValid(mMismatches.getValue(100))) {
              nextStates.add(step);
            } else {
              finalStates.add(step);
            }
          }
        }
      }
      Collections.sort(nextStates, new MismatchComparator());

      final ArrayList<AlignmentState> filteredState  = new ArrayList<>();
      nextStates: for (AlignmentState state : nextStates) {
        for (AlignmentState filtered : filteredState) {
          if (filtered.mChain.mContigId == state.mChain.mContigId
              && filtered.mContigPosition == state.mContigPosition
              && filtered.numberMismatches() < state.numberMismatches()) {
            continue nextStates;
          }
        }
        filteredState.add(state);
      }
      states = filteredState;
    }
    return finalStates;
  }

  /**
   * Entry point
   * @param args command line args
   * @throws IOException should reading fail
   */
  public static void main(String[] args) throws IOException {
    final SequencesReader reader = SequencesReaderFactory.createDefaultSequencesReader(new File(args[0]));
    final Graph graph = GraphReader.read(new StoreDirProxy(new File(args[1])));
    final IntegerOrPercentage percentage;
    if (args.length > 2) {
      percentage = IntegerOrPercentage.valueOf(args[2]);
    } else {
      percentage = IntegerOrPercentage.valueOf(3);
    }
    final ApplyReference apply = new ApplyReference(graph, 18, 18, percentage);
    for (int i = 0; i < reader.numberSequences(); i++) {
      apply.referenceSequence(reader.name(i), reader.read(i), System.out);

    }

  }

}
