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

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import com.rtg.assembler.graph.Contig;
import com.rtg.assembler.graph.Graph;
import com.rtg.assembler.graph.PathsIterator;
import com.rtg.mode.DNA;
import com.rtg.mode.DnaUtils;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.Histogram;
import com.rtg.util.Pair;
import com.rtg.util.StringUtils;
import com.rtg.util.array.intindex.IntChunks;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 */
public class DeBruijnGraphBuilderTest extends TestCase {

  static List<Long> intLinks(long ... links) {
    final List<Long> list = new ArrayList<>();
    for (final long l : links) {
      list.add(l);
    }
    return list;
  }

  //TODO add test for k > 32

  public void testSize() throws IOException {
    Diagnostic.setLogStream();
    final File tmpDir0 =  ReaderTestUtils.getDNADir(">read1" + StringUtils.LS + "accgttc" + StringUtils.LS);
    final File tmpDir1 =  ReaderTestUtils.getDNADir(">read2" + StringUtils.LS + "ccgttcg" + StringUtils.LS
        );
    try {
      final List<ReadPairSource> files = new LinkedList<>();
      files.add(ReadPairSource.makeSource(tmpDir0, LongRange.NONE));
      files.add(ReadPairSource.makeSource(tmpDir1, LongRange.NONE));
      try {
        assertEquals(10, DeBruijnGraphBuilder.size(files, 3));
        assertEquals(2, DeBruijnGraphBuilder.size(files, 7));
        assertEquals(0, DeBruijnGraphBuilder.size(files, 8));
      } finally {
        closeSources(files);
      }
    } finally {
      FileHelper.deleteAll(tmpDir0);
      FileHelper.deleteAll(tmpDir1);
    }
  }

  private void closeSources(List<ReadPairSource> files) {
    for (ReadPairSource r : files) {
      r.close();
    }
  }

  public void testBuild() throws IOException {
    final File tmpDir =  ReaderTestUtils.getDNADir(">read1" + StringUtils.LS + "accgttc" + StringUtils.LS
        + ">read2" + StringUtils.LS + "ccgttcg" + StringUtils.LS
        );
    try {
      final List<ReadPairSource> files = Collections.singletonList(ReadPairSource.makeSource(tmpDir, LongRange.NONE));
      try {
        assertEquals(10, DeBruijnGraphBuilder.size(files, 3));
        assertEquals(2, DeBruijnGraphBuilder.size(files, 7));
        assertEquals(0, DeBruijnGraphBuilder.size(files, 8));
        final DeBruijnGraphBuilder dbg = new DeBruijnGraphBuilder(files, 3, 3, 4);
        final DeBruijnGraph graph = dbg.mDeBruijnGraph;
        graph.setThreshold(0);
        final Map<String, Integer> stringified = new HashMap<>();
        for (final Kmer k: graph) {
          stringified.put(k.toString(), graph.frequency(k));
        }
        assertEquals(1, stringified.get("ACC").intValue());
        assertEquals(2, stringified.get("CCG").intValue());
        assertEquals(null, stringified.get("GCC"));

        assertEquals(2, stringified.get("GAA").intValue());
        assertEquals(null, stringified.get("TTC"));
        // Last kmer of read 2. only once in reverse complement
        assertEquals(1, stringified.get("CGA").intValue());
        assertEquals(null, stringified.get("TCG"));
      } finally {
        closeSources(files);
      }
    } finally {
      FileHelper.deleteAll(tmpDir);
    }
  }
  public void testTwoContigs() throws IOException {
    for (final KmerFactory fact : new KmerFactory[] {StringKmer.factory(), ByteKmer.factory()}) {
      final File tmpDir =  ReaderTestUtils.getDNADir(
          ">read1" + StringUtils.LS + "ataaat" + StringUtils.LS
          + ">read1b" + StringUtils.LS + "ataaat" + StringUtils.LS
          + ">read2" + StringUtils.LS + "ataaac" + StringUtils.LS
          + ">read2b" + StringUtils.LS + "ataaac" + StringUtils.LS
          );
      try {
        try (ReadPairSource o = ReadPairSource.makeSource(tmpDir, LongRange.NONE)) {
          final DeBruijnGraphBuilder dbg = new DeBruijnGraphBuilder(Collections.singletonList(o), 4, fact, 3, 4);
          dbg.buildPreContigs();

          final Map<Long, ComparisonNode> checker = new HashMap<>();
          checker.put(1L, new ComparisonNode("ATAAA", intLinks(), intLinks(2, 3)));
          checker.put(2L, new ComparisonNode("AAAT", intLinks(1), intLinks()));
          checker.put(3L, new ComparisonNode("AAAC", intLinks(1), intLinks()));
          compareGraph(checker, dbg.preContigGraph());
        }

      } finally {
        FileHelper.deleteAll(tmpDir);
      }
    }
  }
  public void testPrecontigs() throws IOException {
    for (final KmerFactory fact : new KmerFactory[] {StringKmer.factory(), ByteKmer.factory()}) {
      // see graph1.jpg
      final File tmpDir =  ReaderTestUtils.getDNADir(
          ">read1" + StringUtils.LS + "ataaat" + StringUtils.LS
          + ">read1b" + StringUtils.LS + "ataaat" + StringUtils.LS
          + ">read2" + StringUtils.LS + "ataaacatccgaat" + StringUtils.LS
          + ">read2b" + StringUtils.LS + "ataaacatccgaat" + StringUtils.LS
          + ">read3" + StringUtils.LS + "ataaacagccgaat" + StringUtils.LS
          + ">read3b" + StringUtils.LS + "ataaacagccgaatc" + StringUtils.LS
          );
      try {
        try (ReadPairSource o = ReadPairSource.makeSource(tmpDir, LongRange.NONE)) {
          final DeBruijnGraphBuilder dbg = new DeBruijnGraphBuilder(Collections.singletonList(o), 4, fact, 0, 4);
          dbg.setGoodThreshold(1);
          dbg.buildPreContigs();

          final Map<Long, ComparisonNode> checker = new HashMap<>();
          checker.put(1L, new ComparisonNode("ATAAA", intLinks(), intLinks(2, 3)));
          checker.put(2L, new ComparisonNode("AAAT", intLinks(1), intLinks()));
          checker.put(3L, new ComparisonNode("AAACA", intLinks(1), intLinks(4, 5)));
          checker.put(4L, new ComparisonNode("ACATCCG", intLinks(3), intLinks(6)));
          checker.put(5L, new ComparisonNode("ACAGCCG", intLinks(3), intLinks(6)));
          checker.put(6L, new ComparisonNode("CCGAAT", intLinks(4, 5), intLinks()));
          compareGraph(checker, dbg.preContigGraph());
          final int[] expectedStartTipValues = {5, 6, 7, 11, 11, 14};
          final int[] expectedEndTipValues = {14, 4, 12, 10, 10, 6};
          checkTipValues(expectedStartTipValues, expectedEndTipValues, checker, dbg);
          dbg.deleteTips(dbg.calculateTipValues());
          final Map<Long, ComparisonNode> noTips = new HashMap<>();
          noTips.put(4L, new ComparisonNode("ACATCCG", intLinks(), intLinks()));
          noTips.put(5L, new ComparisonNode("ACAGCCG", intLinks(), intLinks()));
          compareGraph(noTips, dbg.preContigGraph());
        }
      } finally {
        FileHelper.deleteAll(tmpDir);
      }
    }
  }

  public void checkTipValues(int[] expectedStartTipValues, int[] expectedEndTipValues, Map<Long, ComparisonNode> graph, DeBruijnGraphBuilder dbg) {
    final Pair<IntChunks, IntChunks> tipValues = dbg.calculateTipValues();
    final List<Long> contigIds = new ArrayList<>(graph.keySet());
    Collections.sort(contigIds);
    for (int i = 0; i < graph.size(); ++i) {
      final ComparisonNode node = graph.get(contigIds.get(i));
      final long comparisonId = node.mComparison;
      final IntChunks start = node.mReverseComparison ? tipValues.getB() : tipValues.getA();
      final IntChunks end = node.mReverseComparison ? tipValues.getA() : tipValues.getB();
      assertEquals(expectedStartTipValues[i], start.get(comparisonId));
      assertEquals(expectedEndTipValues[i], end.get(comparisonId));
    }
  }
  public void testPrecontigs2() throws IOException {
    // see graph2.jpg
    final File tmpDir =  ReaderTestUtils.getDNADir(
        ">read1a" + StringUtils.LS + "aaact" + StringUtils.LS
        + ">read1b" + StringUtils.LS + "aaact" + StringUtils.LS
        + ">read2a" + StringUtils.LS + "aaacag" + StringUtils.LS
        + ">read2b" + StringUtils.LS + "aaacag" + StringUtils.LS
        + ">read3a" + StringUtils.LS + "aacacgg" + StringUtils.LS
        + ">read3b" + StringUtils.LS + "aacacgg" + StringUtils.LS
        + ">read4a" + StringUtils.LS + "aacacc" + StringUtils.LS
        + ">read4b" + StringUtils.LS + "aacacc" + StringUtils.LS
        + ">read5a" + StringUtils.LS + "tacgg" + StringUtils.LS
        + ">read5b" + StringUtils.LS + "tacgg" + StringUtils.LS
        );
    try {
      try (ReadPairSource o = ReadPairSource.makeSource(tmpDir, LongRange.NONE)) {
        final DeBruijnGraphBuilder dbg = new DeBruijnGraphBuilder(Collections.singletonList(o), 4, 3, 4);
        dbg.buildPreContigs();
        final int[] expectedStartTipValues = {4, 5, 5, 6, 6, 7, 7, 4, 8};
        final int[] expectedEndTipValues = {8, 4, 7, 4, 6, 5, 4, 5, 4};
        final Map<Long, ComparisonNode> checker = new HashMap<>();
        checker.put(1L, new ComparisonNode("AAAC", intLinks(), intLinks(2, 3)));
        checker.put(2L, new ComparisonNode("AACT", intLinks(1), intLinks()));
        checker.put(3L, new ComparisonNode("AACA", intLinks(1), intLinks(4, 5)));
        checker.put(4L, new ComparisonNode("ACAG", intLinks(3), intLinks()));
        checker.put(5L, new ComparisonNode("ACAC", intLinks(3), intLinks(6, 7)));
        checker.put(6L, new ComparisonNode("CACG", intLinks(5), intLinks(9)));
        checker.put(7L, new ComparisonNode("CACC", intLinks(5), intLinks()));
        checker.put(8L, new ComparisonNode("TACG", intLinks(), intLinks(9)));
        checker.put(9L, new ComparisonNode("ACGG", intLinks(8, 6), intLinks()));
        compareGraph(checker, dbg.preContigGraph());
        checkTipValues(expectedStartTipValues, expectedEndTipValues, checker, dbg);
      }
    } finally {
      FileHelper.deleteAll(tmpDir);
    }

  }
  public void testKsizePalindrome() throws IOException {
    final File tmpDir =  ReaderTestUtils.getDNADir(
        ">read1a" + StringUtils.LS + "cgggaattggcg" + StringUtils.LS
        + ">read1b" + StringUtils.LS + "cgggaattggcg" + StringUtils.LS
        );
    try {
      final ReadPairSource o = ReadPairSource.makeSource(tmpDir, LongRange.NONE);
      final DeBruijnGraphBuilder dbg = new DeBruijnGraphBuilder(Collections.singletonList(o), 4, 3, 4);
      dbg.buildPreContigs();
      final int[] expectedStartTipValues = {7, 7, 8};
      final int[] expectedEndTipValues = {12, 12, 8};
      final Map<Long, ComparisonNode> checker = new HashMap<>();
      checker.put(1L, new ComparisonNode("CGGGAAT", intLinks(), intLinks(+3, -3)));
      checker.put(2L, new ComparisonNode("CGCCAAT", intLinks(), intLinks(+3, -3)));
      checker.put(3L, new ComparisonNode("AATT", intLinks(+1, +2), intLinks(-1, -2)));
      compareGraph(checker, dbg.preContigGraph());
      checkTipValues(expectedStartTipValues, expectedEndTipValues, checker, dbg);
      o.close();
    } finally {
      FileHelper.deleteAll(tmpDir);
    }
  }

  public void testLongerPalindrome() throws IOException {
    final File tmpDir =  ReaderTestUtils.getDNADir(
        ">read1a" + StringUtils.LS + "ataacatcgtacgatgccgc" + StringUtils.LS
        + ">read1b" + StringUtils.LS + "ataacatcgtacgatgccgc" + StringUtils.LS
        );
    try {
      final ReadPairSource o = ReadPairSource.makeSource(tmpDir, LongRange.NONE);
      final DeBruijnGraphBuilder dbg = new DeBruijnGraphBuilder(Collections.singletonList(o), 4, 3, 4);
      dbg.buildPreContigs();
      final int[] expectedStartTipValues = {7, 7, 16};
      final int[] expectedEndTipValues = {20, 20, 16};
      final Map<Long, ComparisonNode> checker = new HashMap<>();
      checker.put(1L, new ComparisonNode("ATAACAT", intLinks(), intLinks(+3, -3)));
      checker.put(2L, new ComparisonNode("GCGGCAT", intLinks(), intLinks(+3, -3)));
      checker.put(3L, new ComparisonNode("CATCGTACGATG", intLinks(+1, +2), intLinks(-1, -2)));
      compareGraph(checker, dbg.preContigGraph());
      checkTipValues(expectedStartTipValues, expectedEndTipValues, checker, dbg);
      o.close();
    } finally {
      FileHelper.deleteAll(tmpDir);
    }
  }

  public void testHairpin() throws IOException {
    final String sequence = "ataacatcgatgccgc";
    final File tmpDir =  ReaderTestUtils.getDNADir(
        ">read1a" + StringUtils.LS + sequence + StringUtils.LS
        + ">read1b" + StringUtils.LS + sequence + StringUtils.LS
        + ">read1b" + StringUtils.LS + sequence + StringUtils.LS
        );
    try {
      final ReadPairSource o = ReadPairSource.makeSource(tmpDir, LongRange.NONE);
      final DeBruijnGraphBuilder dbg = new DeBruijnGraphBuilder(Collections.singletonList(o), 4, 3, 4);
      dbg.buildPreContigs();
      final int[] expectedStartTipValues = {7 , 7 , 12 };
      final int[] expectedEndTipValues = {16 , 16 , 12 };
      final Map<Long, ComparisonNode> checker = new HashMap<>();
      checker.put(1L, new ComparisonNode("ATAACAT", intLinks(), intLinks(3, -3)));
      checker.put(2L, new ComparisonNode("GCGGCAT", intLinks(), intLinks(3, -3)));
      checker.put(3L, new ComparisonNode("CATCGATG", intLinks(1, 2), intLinks(-1, -2)));
      compareGraph(checker, dbg.preContigGraph());
      checkTipValues(expectedStartTipValues, expectedEndTipValues, checker, dbg);
      o.close();
    } finally {
      FileHelper.deleteAll(tmpDir);
    }
  }

  public void testPalindrome() throws IOException {
    for (final KmerFactory fact : new KmerFactory[] {StringKmer.factory(), ByteKmer.factory()}) {
      final File tmpDir =  ReaderTestUtils.getDNADir(
          ">read1a" + StringUtils.LS + "CGGGAATT" + StringUtils.LS
          + ">read1b" + StringUtils.LS + "CGGGAATT" + StringUtils.LS
          + ">read1b" + StringUtils.LS + "CAATT" + StringUtils.LS
          + ">read1b" + StringUtils.LS + "CAATT" + StringUtils.LS
          );
      try {
        final ReadPairSource o = ReadPairSource.makeSource(tmpDir, LongRange.NONE);
        final DeBruijnGraphBuilder dbg = new DeBruijnGraphBuilder(Collections.singletonList(o), 4, fact, 3, 4);
        dbg.buildPreContigs();
        final int[] expectedStartTipValues = {7, 4, 8};
        final int[] expectedEndTipValues = {12, 9, 8};
        final Map<Long, ComparisonNode> checker = new HashMap<>();
        checker.put(1L, new ComparisonNode("CGGGAAT", intLinks(), intLinks(3, -3)));
        checker.put(2L, new ComparisonNode("CAAT", intLinks(), intLinks(3, -3)));
        checker.put(3L, new ComparisonNode("AATT", intLinks(1, 2), intLinks(-1, -2)));
        compareGraph(checker, dbg.preContigGraph());
        checkTipValues(expectedStartTipValues, expectedEndTipValues, checker, dbg);
        o.close();
      } finally {
        FileHelper.deleteAll(tmpDir);
      }
    }
  }
  public void testDoublePalindrome() throws IOException {
    for (final KmerFactory fact : new KmerFactory[] {StringKmer.factory(), ByteKmer.factory()}) {
      final String sequence = "ACGAAATTGGCCCTAC";
      final File tmpDir =  ReaderTestUtils.getDNADir(
          ">read1a" + StringUtils.LS + sequence + StringUtils.LS
          + ">read1b" + StringUtils.LS + sequence + StringUtils.LS
          );
      try {
        final ReadPairSource o = ReadPairSource.makeSource(tmpDir, LongRange.NONE);
        final DeBruijnGraphBuilder dbg = new DeBruijnGraphBuilder(Collections.singletonList(o), 4, fact, 3, 4);
        dbg.buildPreContigs();
        final int[] expectedStartTipValues = {7, 0, 0, 0, 0};
        final int[] expectedEndTipValues = {0, 7, 0, 0, 0};
        final Map<Long, ComparisonNode> checker = new HashMap<>();
        checker.put(1L, new ComparisonNode("ACGAAAT", intLinks(), intLinks(3, -3)));
        checker.put(2L, new ComparisonNode("GCCCTAC", intLinks(4, -4), intLinks()));
        checker.put(3L, new ComparisonNode("AATT", intLinks(1, -5), intLinks(5, -1)));
        checker.put(4L, new ComparisonNode("GGCC", intLinks(5, -2), intLinks(-5, 2)));
        checker.put(5L, new ComparisonNode("ATTGGC", intLinks(3, -3), intLinks(4, -4)));
        compareGraph(checker, dbg.preContigGraph());
        checkTipValues(expectedStartTipValues, expectedEndTipValues, checker, dbg);
        o.close();
      } finally {
        FileHelper.deleteAll(tmpDir);
      }
    }
  }

  public void testDoublePalindromeOuter() throws IOException {
    // Slightly dodgy test. Results are dependant on the hashing order
    for (final KmerFactory fact : new KmerFactory[] {StringKmer.factory()}) {
      final String sequence = "AATTGCTCGGCC";
      final File tmpDir =  ReaderTestUtils.getDNADir(
          ">read1a" + StringUtils.LS + sequence + StringUtils.LS
          + ">read1b" + StringUtils.LS + sequence + StringUtils.LS
          );
      try {
        final ReadPairSource o = ReadPairSource.makeSource(tmpDir, LongRange.NONE);
        final DeBruijnGraphBuilder dbg = new DeBruijnGraphBuilder(Collections.singletonList(o), 4, fact, 3, 4);
        dbg.buildPreContigs();
        assertEquals(1, dbg.preContigGraph().numberContigs());
        assertEquals(19, dbg.preContigGraph().contigLength(1));
        o.close();
      } finally {
        FileHelper.deleteAll(tmpDir);
      }
    }
  }

  static String contigAsString(Contig c) {
    final StringBuilder sb = new StringBuilder();
    for (int i = 0; i < c.length(); ++i) {
      sb.append(DNA.valueChars()[c.nt(i)]);
    }
    return sb.toString();
  }
  static void compareGraph(Map<Long, ComparisonNode> expected, Graph graph) {
    //    assertEquals(expected.size(), graph.numberContigs());
    int foundCount = 0;
    for (long i = 1; i <= graph.numberContigs(); ++i) {
      if (graph.contigDeleted(i)) {
        continue;
      }
      final Contig contig = graph.contig(i);
      final String sequence = contigAsString(contig);
      boolean found = false;
      for (final ComparisonNode expectedNode : expected.values()) {
        if (expectedNode.matches(sequence)) {
          if (expectedNode.getComparison() != -1) {
            fail("Node appears twice:" + sequence + " " + i + " & " + expectedNode.mComparison);
          }
          expectedNode.setComparison(i, sequence);
          found = true;
          ++foundCount;
        }
      }
      assertTrue("We weren't expecting to see: " + sequence, found);
    }
    assertEquals(expected.size(), foundCount);
    for (final ComparisonNode expectedNode : expected.values()) {
      final Pair<List<Long>, List<Long>> links = actualLinks(expectedNode.mComparison * (expectedNode.mReverseComparison ? -1 : 1), graph);
      assertEquals(actualIds(expectedNode.mStartLinks, expected), links.getA());
      assertEquals(actualIds(expectedNode.mEndLinks, expected), links.getB());
    }

  }
  static List<Long> actualIds(List<Long> expectedIds, Map<Long, ComparisonNode> map) {
    final List<Long> exp = new ArrayList<>();
    for (final long l : expectedIds) {
      final ComparisonNode comparisonNode = map.get(Math.abs(l));
      exp.add(comparisonNode.mComparison * (comparisonNode.mReverseComparison ? -1 : 1) * (l < 0 ? -1 : 1));
    }
    Collections.sort(exp);
    return exp;
  }
  static Pair<List<Long>, List<Long>> actualLinks(long contigId, Graph actual) {
    final List<Long> startLinks = new ArrayList<>();
    final List<Long> endLinks = new ArrayList<>();
    final PathsIterator actualPaths = actual.paths(contigId);
    long pathId;
    while ((pathId = actualPaths.nextPathId()) != 0) {
      assert actual.pathLength(pathId) == 2;

      final int index = actualPaths.contigIndex();
      final long linkedContig = actual.pathContig(pathId, 1 - index);
      if (index == 0) {
        endLinks.add(linkedContig);
      } else {
        startLinks.add(linkedContig);
      }
    }
    Collections.sort(startLinks);
    Collections.sort(endLinks);
    return Pair.create(startLinks, endLinks);
  }

  static class ComparisonNode {
    final String mSequence;
    final List<Long> mStartLinks;
    final List<Long> mEndLinks;
    private long mComparison = -1;
    boolean mReverseComparison;

    ComparisonNode(String sequence, List<Long> startLinks, List<Long> endLinks) {
      mSequence = sequence;
      mStartLinks = startLinks;
      mEndLinks = endLinks;
    }

    boolean matches(String nodeSequence) {
      return nodeSequence.equals(mSequence) || nodeSequence.equals(DnaUtils.reverseComplement(mSequence));
    }
    long getComparison() {
      return mComparison;
    }
    void setComparison(long id, String sequence) {
      mComparison = id;
      mReverseComparison = !sequence.equals(mSequence);
    }

    @Override
    public String toString() {
      return "ComparisonNode{"
          + "mSequence='" + mSequence + '\''
          + ", mStartLinks=" + mStartLinks
          + ", mEndLinks=" + mEndLinks
          + ", mComparison=" + mComparison
          + ", mReverseComparison=" + mReverseComparison
          + '}';
    }
  }

  // From Rhodobacter sphaeroides  D0F7EACXX
  //private static final int[] HD0F7 = {0, 63746955, 2437194, 641039, 559275, 556284, 525495, 466861, 396736, 327456, 267143, 221695, 191280, 168788, 154597, 145237, 138523, 133691, 129294, 125222, 121238, 117619, 115752, 112053, 107534, 102632, 98907, 94030, 88786, 85080, 78776, 73679, 69656, 65153, 60582, 55751, 51709, 46995, 42353, 38171, 34825, 31462, 27713, 24888, 21515, 19183, 16702, 14718, 12744, 11207, 9786, 8754, 7561, 6764, 5960, 5256, 4691, 4106, 3802, 3376, 3271, 2965, 2648, 2419, 2206, 2027, 1848, 1813, 1715, 1676, 1466, 1464, 1346, 1306, 1165, 1166, 1174, 1147, 1030, 977, 957, 896, 838, 791, 819, 749, 747, 715, 715, 607, 606, 574, 554, 486, 489, 504, 454, 465, 445, 434, 422, 404, 361, 368, 343, 324, 317, 288, 262, 303, 230, 218, 211, 229, 201, 204, 200, 205, 177, 156, 173, 139, 148, 117, 134, 121, 124, 119, 128, 107, 121, 110, 117, 103, 101, 103, 103, 81, 78, 95, 102, 90, 62, 86, 91, 70, 64, 63, 60, 57, 56, 64, 73, 46, 49, 38, 31, 33, 31, 32, 23, 31, 24, 16, 22, 15, 19, 17};

  // From Rhodobacter sphaeroides gerald_81RE8ABXX_3_GCCAAT
  private static final int[] H81RE = {0, 342341000, 21135180, 7829310, 3807240, 2056229, 1190860, 729557, 468777, 314810, 219507, 159911, 120109, 92156, 72979, 58405, 48534, 40519, 34093, 29384, 25409, 22735, 20235, 17933, 16403, 15250, 14113, 13706, 12792, 12271, 11917, 11477, 11239, 10988, 10489, 10348, 10288, 10499, 10165, 10123, 10325, 10408, 10490, 10482, 10811, 10676, 10924, 11042, 10882, 11346, 11113, 11345, 11527, 11811, 11995, 12626, 12718, 12882, 13172, 13290, 13297, 13489, 13677, 13883, 13982, 14123, 14614, 14543, 14848, 15110, 15208, 15622, 15788, 16155, 16403, 16541, 16967, 16863, 17630, 17845, 17883, 17980, 18297, 18415, 18863, 18810, 19116, 19338, 19386, 19980, 19986, 19914, 20249, 20477, 20697, 21105, 21115, 21818, 21533, 21808, 21833, 22027, 22553, 22764, 22863, 23059, 23274, 23685, 23656, 23867, 24126, 24280, 24218, 24695, 24394, 24628, 24675, 24991, 24833, 24874, 24967, 25689, 25553, 25450, 25792, 25720, 25848, 25759, 25515, 25627, 25954, 26044, 26382, 26160, 26460, 25987, 26240, 26297, 26616, 26642, 26490, 26547, 26496, 26450, 26301, 26052, 26567, 26684, 26470, 26000, 26004, 25969, 26037, 25880, 25787, 25412, 25674, 25339, 25296, 25219, 25485, 25175, 24855, 25128, 24856, 24719, 24449, 24642, 24398, 23987, 23823, 23764, 23637, 23141, 22973, 22671, 23003, 22699, 22525, 22120, 21956, 21902, 21705, 21724, 21119, 21120, 20855, 20614, 20256, 19964, 19980, 19710, 19514, 19346, 19089, 18764, 18633, 18660, 18273, 18460, 17809, 17942, 17827, 17410, 17126, 17082, 16803, 16524, 15986, 15989, 15774, 15701, 15071, 14905, 14756, 14684, 14295, 14208, 13911, 13804, 13712, 13260, 13189, 12878, 12769, 12555, 12397, 12070, 11845, 11841, 11546, 11488, 11056, 10989, 10816, 10596, 10399, 10214, 10115, 9880, 9681, 9420, 9248, 9075, 8914, 8730, 8453, 8271, 8195, 8005, 8001, 7624, 7502, 7332, 7244, 7197, 6900, 6921, 6825, 6443, 6494, 6248, 6115, 6136, 5875, 5841, 5679, 5411, 5265, 5195, 5187, 4959, 4983, 4793, 4742, 4674, 4561, 4443, 4323, 4284, 4123, 4011, 3952, 3819, 3800, 3590, 3640, 3469, 3324, 3230, 3223, 3093, 3014, 2879, 2795, 2750, 2719, 2716, 2668, 2572, 2545, 2479, 2457, 2395, 2273, 2310, 2313, 2225, 2203, 2083, 2030, 1996, 2002, 1882, 1973, 1754, 1849, 1780, 1715, 1657, 1617, 1694, 1590, 1607, 1419, 1411, 1381, 1379, 1321, 1297, 1266, 1201, 1186, 1111, 1127, 1111, 1103, 1037, 931, 942, 955, 924, 905, 871, 887, 872, 873, 860, 823, 819, 800, 757, 749, 725, 719, 720, 709, 701, 729, 723, 685, 671, 634, 614, 615, 575, 624, 582, 575, 548, 546, 491, 504, 504, 454, 460, 521, 476, 408, 478, 449, 391, 420, 422, 422, 402, 381, 381, 386, 394, 368, 381, 429, 383, 396, 364, 379, 362, 356, 397, 368, 360, 357, 348, 329, 365, 362, 347, 342, 309, 335, 336, 306, 347, 282, 326, 337, 307, 325, 256, 293, 263, 302, 276, 275, 288, 262, 276, 250, 257, 246, 294, 282, 234, 257, 279, 252, 257, 212, 260, 263, 214, 217, 236, 254, 257, 241, 245, 239, 237, 232, 262, 258, 277, 252, 236, 253, 254, 243, 233, 246, 258, 228, 257, 244, 240, 236, 226, 248, 206, 252, 196, 213, 228, 222, 208, 213, 220, 212, 199, 195, 198, 217, 208, 216, 223, 195, 211, 194, 223, 181, 197, 182, 191, 223, 197, 206, 183, 192, 204, 176, 198, 222, 200, 190, 211, 171, 203, 186, 154, 178, 189, 171, 220, 181, 166, 161, 168, 180, 159, 174, 182, 216, 175, 176, 181, 168, 181, 180, 193, 165, 187, 185, 176, 166, 147, 177, 175, 187, 160, 182, 161, 171, 176, 162, 157, 156, 187, 159, 156, 143, 171, 173, 129, 153, 139, 170, 150, 147, 141, 126, 148, 129, 152, 124, 130, 158, 126, 119, 115, 124, 129, 129, 118, 116, 101, 107, 114, 122, 98, 115, 111, 109, 111, 82, 97, 100, 103, 103, 114, 90, 87, 94, 86, 95, 94, 76, 90, 81, 100, 89, 90, 92, 84, 82, 89, 71, 84, 76, 83, 79, 85, 78, 79, 78, 67, 76, 69, 71, 72, 93, 92, 69, 78, 69, 77, 67, 73, 73, 64, 64, 59, 63, 75, 66, 64, 62, 67, 66, 66, 67, 62, 69, 67, 69, 51, 86, 77, 57, 69, 83, 62, 68, 60, 69, 60, 49, 75, 56, 69, 84, 80, 54, 60, 67, 60, 62, 60, 62, 65, 74, 77, 61, 68, 72, 53, 66, 71, 55, 68, 51, 62, 68, 52, 57, 63, 58, 43, 43, 55, 71, 62, 40, 52, 63, 58, 60, 47, 38, 55, 46, 59, 50, 43, 51, 50, 55, 58, 44, 50, 44, 51, 43, 39, 48, 41, 41, 53, 56, 43, 39, 39, 33, 42, 45, 45, 52, 40, 40, 46, 47, 35, 40, 43, 42, 32, 38, 48, 45, 38, 29, 45, 33, 29, 40, 38, 28, 35, 41, 29, 38, 35, 36, 29, 36, 43, 35, 35, 34, 35, 27, 28, 28, 30, 35, 37, 36, 37, 27, 30, 25, 33, 27, 28, 39, 27, 31, 30, 30, 30, 30, 42, 33, 36, 38, 36, 33, 37, 31, 27, 30, 47, 35, 35, 27, 34, 40, 39, 33, 28, 30, 38, 29, 42, 39, 27, 33, 26, 28, 27, 24, 21, 41, 18, 31, 27, 35, 25, 29, 26, 37, 26, 29, 33, 23, 20, 27, 20, 18, 28, 35, 31, 22, 24, 26, 28, 34, 16, 43, 23, 27, 29, 18, 24, 33, 24, 19, 23, 32, 19, 21, 32, 28, 33, 21, 27, 23, 26, 26, 21, 30, 33, 16, 31, 21, 24, 33, 25, 24, 34, 28, 32, 28, 38, 25, 29, 23, 28, 29, 30, 30, 23, 21, 24, 26, 18, 27, 20, 22, 27, 31, 26, 21, 19, 15, 20, 23, 12, 16, 18, 20, 16, 10, 19, 24, 21, 20, 11, 31, 18, 13, 15, 22, 17, 18, 24, 15, 18, 16, 15, 20, 21, 19, 13, 15, 21, 17, 15, 12, 20, 15, 9, 15, 9, 13, 13, 13, 14, 20, 10, 11, 16, 16, 9, 11, 14, 13, 12, 13, 14, 19, 14, 12, 15, 12, 12, 15, 10, 11, 12, 20, 14, 16, 15, 13, 11, 12, 16, 12, 6, 17, 15, 17, 23, 12, 7, 16, 16, 16, 19, 11, 26, 14, 15, 19, 17, 19, 13, 16, 16, 20, 16, 10, 14, 12, 4, 14, 9, 12, 8, 2744};

  // A real Enterococcus faecalis distribution
  private static final int[] HEF = {0, 68649596, 4825528, 1011684, 359389, 176421, 102581, 67024, 45800, 33434, 24832, 19090, 14745, 12122, 10007, 8077, 6877, 5656, 4818, 4147, 3690, 3217, 2900, 2461, 2184, 1965, 1782, 1601, 1420, 1286, 1270, 1073, 1023, 901, 895, 777, 749, 646, 662, 565, 517, 500, 459, 449, 404, 338, 320, 326, 304, 281, 255, 298, 264, 242, 208, 172, 237, 188, 166, 158, 175, 167, 163, 152, 133, 130, 130, 124, 102, 103, 97, 114, 119, 90, 95, 91, 88, 88, 107, 85, 92, 73, 75, 81, 61, 66, 63, 73, 71, 62, 60, 61, 76, 57, 57, 62, 53, 54, 56, 52, 66, 53, 59, 55, 47, 51, 50, 65, 64, 57, 58, 69, 62, 60, 66, 70, 88, 85, 69, 93, 67, 78, 82, 90, 86, 106, 92, 99, 106, 116, 117, 103, 115, 107, 113, 104, 101, 106, 108, 108, 115, 108, 132, 134, 135, 154, 121, 143, 148, 145, 177, 157, 198, 211, 229, 220, 219, 209, 215, 226, 252, 264, 290, 267, 287, 292, 310, 350, 323, 357, 386, 351, 393, 426, 420, 407, 441, 452, 460, 482, 520, 529, 533, 544, 581, 613, 628, 643, 631, 684, 681, 685, 701, 774, 762, 785, 819, 838, 881, 941, 954, 1036, 994, 1022, 1053, 1209, 1236, 1265, 1342, 1319, 1432, 1538, 1598, 1635, 1697, 1750, 1787, 1881, 1951, 2066, 2219, 2191, 2379, 2432, 2526, 2714, 2755, 2868, 2978, 3129, 3114, 3239, 3527, 3543, 3713, 3946, 3956, 4254, 4442, 4628, 4801, 4844, 5256, 5275, 5517, 5635, 5971, 6143, 6221, 6564, 6699, 6950, 7517, 7780, 7846, 8249, 8613, 8804, 9131, 9344, 9666, 9911, 10130, 10610, 10736, 11251, 11491, 11882, 12253, 12452, 12830, 13301, 13705, 13955, 14045, 14627, 14863, 15138, 15636, 15547, 15727, 16381, 16466, 17006, 16968, 17592, 17802, 18299, 18492, 18912, 19062, 19432, 19588, 19767, 19953, 20543, 20480, 20460, 20860, 21205, 21263, 21498, 21381, 21923, 22354, 22453, 22391, 22614, 22147, 23038, 22457, 23060, 22971, 23198, 23259, 23146, 23245, 23311, 23174, 23374, 22912, 23254, 23202, 23284, 23307, 22987, 22992, 23273, 23051, 22919, 22919, 22705, 22801, 22626, 22269, 22243, 22208, 22382, 21938, 21843, 21929, 21703, 21884, 21494, 21425, 21173, 21087, 21492, 21124, 21031, 20941, 20585, 20721, 20122, 20076, 20027, 19816, 19559, 19570, 19401, 19041, 18994, 18547, 18415, 18245, 17966, 17797, 17590, 17176, 16803, 16887, 16643, 16330, 16192, 15834, 15813, 15580, 15426, 15052, 14884, 14610, 14358, 14178, 13936, 13749, 13298, 13251, 12581, 12552, 12323, 12074, 11602, 11398, 11163, 10925, 10772, 10517, 10208, 9978, 9844, 9445, 9038, 8538, 8586, 8450, 8204, 7945, 7711, 7334, 7287, 7004, 6959, 6676, 6428, 6326, 6031, 6019, 5784, 5538, 5219, 5003, 4954, 4781, 4675, 4569, 4199, 4125, 3920, 3846, 3625, 3519, 3401, 3198, 3263, 3088, 3054, 2836, 2791, 2719, 2492, 2444, 2405, 2254, 2142, 2040, 1937, 1875, 1770, 1731, 1571, 1561, 1501, 1433, 1335, 1393, 1302, 1170, 1187, 1089, 1084, 1025, 931, 920, 870, 815, 792, 685, 689, 675, 619, 534, 513, 481, 485, 490, 445, 414, 403, 383, 401, 337, 337, 315, 327, 293, 267, 223, 235, 211, 220, 256, 242, 186, 178, 167, 184, 180, 165, 153, 152, 139, 125, 139, 113, 134, 125, 130, 115, 102, 91, 92, 90, 83, 89, 92, 85, 85, 62, 65, 67, 35, 50, 54, 51, 54, 53, 49, 62, 48, 46, 45, 54, 53, 43, 26, 38, 52, 39, 44, 40, 47, 39, 42, 42, 32, 38, 32, 38, 32, 35, 33, 32, 27, 38, 38, 33, 43, 36, 33, 37, 28, 30, 35, 31, 34, 35, 37, 40, 39, 33, 40, 27, 43, 50, 39, 30, 31, 42, 36, 43, 47, 40, 51, 37, 44, 51, 49, 46, 46, 41, 39, 56, 43, 53, 55, 41, 60, 50, 55, 53, 51, 54, 56, 56, 49, 48, 58, 67, 61, 87, 75, 72, 72, 82, 64, 76, 88, 71, 78, 83, 83, 77, 90, 66, 78, 83, 75, 84, 97, 80, 102, 77, 92, 80, 83, 95, 84, 92, 101, 94, 86, 109, 88, 87, 114, 113, 114, 101, 123, 101, 94, 111, 116, 134, 116, 119, 102, 126, 131, 134, 132, 119, 142, 122, 159, 152, 120, 132, 139, 122, 125, 152, 141, 126, 122, 127, 128, 116, 139, 134, 141, 167, 161, 159, 163, 180, 164, 160, 156, 182, 183, 158, 176, 175, 185, 185, 167, 168, 176, 178, 195, 165, 170, 171, 177, 221, 190, 191, 187, 210, 176, 194, 183, 186, 177, 192, 150, 180, 175, 194, 168, 196, 189, 167, 184, 181, 180, 173, 173, 187, 195, 182, 199, 208, 175, 177, 207, 185, 183, 190, 212, 229, 229, 214, 215, 212, 204, 220, 241, 206, 215, 195, 194, 196, 199, 194, 209, 183, 208, 200, 180, 177, 202, 168, 164, 178, 191, 178, 158, 181, 177, 145, 155, 160, 179, 158, 174, 191, 158, 145, 171, 146, 140, 151, 164, 127, 143, 167, 145, 148, 138, 140, 120, 142, 135, 118, 121, 151, 116, 129, 108, 140, 134, 107, 137, 118, 132, 132, 108, 106, 128, 123, 141, 122, 109, 113, 115, 112, 137, 130, 112, 118, 109, 103, 126, 108, 96, 102, 102, 113, 98, 116, 109, 114, 106, 83, 71, 70, 95, 95, 84, 70, 82, 69, 77, 59, 62, 79, 58, 59, 66, 52, 47, 55, 51, 53, 36, 59, 43, 42, 48, 36, 30, 40, 32, 34, 33, 32, 44, 32, 34, 29, 37, 30, 38, 44, 27, 44, 29, 26, 36, 41, 25, 27, 37, 20, 21, 34, 30, 27, 22, 27, 16, 23, 19, 16, 18, 16, 20, 21, 22, 29, 24, 18, 25, 17, 18, 19, 14, 21, 19, 15, 13, 14, 17, 22, 20, 14, 20, 23, 17, 16, 12, 15, 18, 16, 18, 15, 18, 15, 10, 15, 17, 9, 17, 24, 17, 18, 10, 23, 13, 27, 17, 17, 18, 17, 9, 10, 19, 11, 14, 14, 13, 11, 21, 13, 19, 13, 16, 9, 12, 20, 15, 15, 14, 13, 20, 19, 9, 19, 21, 14, 14, 14, 17, 14, 12, 9, 19, 8, 10, 11, 10, 14, 9, 10, 11, 9, 9, 10, 12, 9, 11, 9, 12, 7, 10, 13, 16, 9, 7098};

  private int check(final int[] hist) {
    final Histogram h = new Histogram();
    for (int k = 0; k < hist.length; ++k) {
      h.increment(k, hist[k]);
    }
    return DeBruijnGraphBuilder.computeThreshold(h);
  }

  /*
  public void testThresholdLow() {
    assertEquals(1, check(HD0F7));
  }
  */

  public void testThresholdHigh() {
    assertEquals(18, check(H81RE));
  }

  public void testThresholdEfaecalis() {
    assertEquals(18, check(HEF));
  }
}
