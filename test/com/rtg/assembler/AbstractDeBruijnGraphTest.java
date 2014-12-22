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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import com.rtg.mode.DnaUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.iterators.Transform;

import junit.framework.TestCase;

/**
 */
public abstract class AbstractDeBruijnGraphTest extends TestCase {

  protected static class KmerIterableMock implements KmerIterable {
    private final String[] mKmers;
    KmerIterableMock(final String... kmers) {
      mKmers = kmers;
    }
    private static final Transform<String, Kmer> STRING2KMER = new String2Kmer();

    @Override
    public void close() throws IOException {
    }

    private static class String2Kmer extends Transform<String, Kmer> {
      @Override
      public Kmer trans(String x) {
        return new StringKmer(x);
      }
    }

    @Override
    public Iterator<Kmer> iterator() {
      final Iterator<String> strings = Transform.array2Iterator(mKmers);
      return STRING2KMER.trans(strings);
    }
  }

  protected static class KmerMockFactory implements KmerIterableFactoryInterface {
    private final String[] mKmers;
    KmerMockFactory(final String... kmers) {
      mKmers = kmers;
    }
    @Override
    public KmerIterable makeIterable() {
      return new KmerIterableMock(mKmers);
    }
  }


  /**
   * Construct a string that is bigger than target by concantenating small.
   * @param small small fragment to be multiplied
   * @return a string whose length &ge; target.
   */
  protected static String big(final String small, final int target) {
    final StringBuilder sb = new StringBuilder();
    while (sb.length() < target) {
      sb.append(small);
    }
    return sb.toString();
  }


  abstract DeBruijnGraph getDeBruijnGraph(final KmerMockFactory kit, long size, int kMerSize);

  protected final void aTest(int size, int kmerSize, final String sNot, final String s0, final String s1) {
    Diagnostic.setLogStream();
    assertEquals(kmerSize, sNot.length());
    assertEquals(kmerSize, s0.length());
    assertEquals(kmerSize, s1.length());
    final String s0r = DnaUtils.reverseComplement(s0);
    final String s1r = DnaUtils.reverseComplement(s1);

    final StringKmer notInGraph = new StringKmer(sNot);
    final KmerMockFactory kit = new KmerMockFactory(s0, s0, s0r, s1);
    final DeBruijnGraph graph = getDeBruijnGraph(kit, size, kmerSize);
    checkIterator(graph, s0, s1);
    assertFalse(graph.contains(notInGraph));
    graph.setThreshold(1);
    assertFalse(graph.contains(notInGraph));
    checkIterator(graph, s0);
    assertFalse(graph.contains(notInGraph));
    assertFalse(graph.contains(new StringKmer(s1)));
    graph.setThreshold(3);
    checkIterator(graph);
    assertFalse(graph.contains(notInGraph));
    assertFalse(graph.contains(new StringKmer(s1)));
    assertFalse(graph.contains(new StringKmer(s0)));
    graph.setThreshold(0);
    checkIterator(graph, s0, s1);

    assertEquals(3, graph.frequency(new StringKmer(s0)));
    assertEquals(3, graph.frequency(new StringKmer(s0r)));

    assertEquals(1, graph.frequency(new StringKmer(s1)));
    assertEquals(1, graph.frequency(new StringKmer(s1r)));

    assertFalse(graph.isBuilt(new StringKmer(s0)));
    graph.setBuilt(new StringKmer(s0), true);
    assertTrue(graph.isBuilt(new StringKmer(s0)));
    graph.setBuilt(new StringKmer(s0), false);
    assertFalse(graph.isBuilt(new StringKmer(s0)));

    assertTrue(graph.contains(new StringKmer(s0)));
    assertTrue(graph.contains(new StringKmer(s1r)));
    assertFalse(graph.contains(notInGraph));

    final HashSet<String> expected = new HashSet<>();
    expected.add(s0);
    expected.add(s1);
    final HashSet<String> actual = new HashSet<>();
    for (final Kmer k : graph) {
      actual.add(k.toString());
    }
    assertEquals(expected, actual);

    expected.add(s1);
    actual.clear();
    graph.setThreshold(0);

    assertTrue(graph.contains(new StringKmer(s0)));
    assertTrue(graph.contains(new StringKmer(s0r)));
    assertTrue(graph.contains(new StringKmer(s1r)));
    assertTrue(graph.contains(new StringKmer(s1)));
    assertFalse(graph.contains(notInGraph));
    checkIterator(graph, s0, s1);

    for (final Kmer k : graph) {
      actual.add(k.toString());
    }
    assertEquals(expected, actual);

    try {
      assertEquals(0, graph.frequency(notInGraph));
      fail();
    } catch (final RuntimeException e) {
      //expected
    }
    try {
      assertFalse(graph.isBuilt(notInGraph));
      fail();
    } catch (final RuntimeException e) {
      //expected
    }
    try {
      graph.setBuilt(notInGraph, true);
      fail();
    } catch (final RuntimeException e) {
      //expected
    }
  }

  private void checkIterator(DeBruijnGraph graph, String... strings) {
    final Set<String> exp = new HashSet<>();
    Collections.addAll(exp, strings);
    final Set<String> actual = new HashSet<>();
    final Iterator<Kmer> it = graph.iterator();
    while (it.hasNext()) {
      final Kmer aGraph = it.next();
      assertTrue(graph.contains(aGraph));
      actual.add(aGraph.toString());
    }
    assertFalse(it.hasNext());
    assertEquals(exp, actual);

    for (final String str : strings) {
      assertTrue(graph.contains(new StringKmer(str)));
    }
  }

  public void testAllAs() {
    final DeBruijnGraph deBruijnGraph = getDeBruijnGraph(new KmerMockFactory("TTTTT", "AAAAA", "ACCCA", "AAAAA", "AAAAA", "ACCCC", "ACCCC", "ACCCA", "GGGGT"), 9, 5);
    final List<String> kmers = new ArrayList<>();
    for (Kmer k : deBruijnGraph) {
      kmers.add(k.toString());
    }
    Collections.sort(kmers);
    assertEquals(Arrays.asList("AAAAA", "ACCCA", "ACCCC"), kmers);
  }
}
