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

package com.rtg.assembler.graph.implementation;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

import com.rtg.assembler.graph.Contig;
import com.rtg.assembler.graph.Graph;
import com.rtg.assembler.graph.MutableGraph;
import com.rtg.assembler.graph.Path;
import com.rtg.assembler.graph.PathsIterator;
import com.rtg.util.integrity.Exam;
import com.rtg.variant.dna.DNARange;

import junit.framework.TestCase;


/**
 */
public class GraphImplementationTest extends TestCase {

  //normal case for adds then test attribute setting and getting
  public void testAttributes() {
    final GraphImplementation graph = graph();

    assertEquals(null, graph.contigAttribute(1, "foo"));
    assertEquals(null, graph.contigAttribute(-1, "foo"));
    assertEquals(null, graph.contigAttribute(1, "bar"));

    graph.setContigAttribute(1, "foo", "baz");
    // We are going to delete contig 4 before compacting so give contig 5 an attribute
    graph.setContigAttribute(5, "foo", "bang");
    assertEquals("baz", graph.contigAttribute(1, "foo"));
    assertEquals("baz", graph.contigAttribute(-1, "foo"));
    assertEquals(null, graph.contigAttribute(2, "foo"));
    assertEquals(null, graph.contigAttribute(1, "bar"));
    assertEquals("bang", graph.contigAttribute(5, "foo"));

    try {
      graph.contigAttribute(1, "xxx");
      fail();
    } catch (final IllegalArgumentException e) {
      // expected
    }
    try {
      graph.contigAttribute(0, "foo");
      fail();
    } catch (final IllegalArgumentException e) {
      // expected
    }

    assertEquals(null, graph.pathAttribute(1, "boo"));
    assertEquals(null, graph.pathAttribute(-1, "boo"));
    assertEquals(null, graph.pathAttribute(1, "oob"));

    graph.setPathAttribute(1, "boo", "zab");
    assertEquals("zab", graph.pathAttribute(1, "boo"));
    assertEquals("zab", graph.pathAttribute(-1, "boo"));
    assertEquals(null, graph.pathAttribute(2, "boo"));
    assertEquals(null, graph.pathAttribute(1, "oob"));

    try {
      graph.pathAttribute(1, "xxx");
      fail();
    } catch (final IllegalArgumentException e) {
      // expected
    }
    try {
      graph.pathAttribute(0, "foo");
      fail();
    } catch (final IllegalArgumentException e) {
      // expected
    }

    //graph.dump(System.err);
    graph.deleteContig(4);
    final MutableGraph compact = graph.compact();
    assertEquals("baz", compact.contigAttribute(1, "foo"));
    assertEquals("baz", compact.contigAttribute(-1, "foo"));
    assertEquals(null, compact.contigAttribute(2, "foo"));
    assertEquals(null, compact.contigAttribute(1, "bar"));

    assertEquals("zab", compact.pathAttribute(1, "boo"));
    assertEquals("zab", compact.pathAttribute(-1, "boo"));
    // Contig 5's attribute has been shifted to contig 4 (id's are displaced due to deletion)
    assertEquals("bang", compact.contigAttribute(4, "foo"));
    assertEquals(null, compact.pathAttribute(2, "boo"));
    assertEquals(null, compact.pathAttribute(1, "oob"));
  }

  //normal case for adds then compact at end
  public void testCompact() {
    final String[] contigNt = {"ACGT", "CGTA", "GTAC", "TACG", "AA", "CC", "GG", "TT"};
    final Contig[] contigs = new Contig[contigNt.length];
    for (int i = 0; i < contigNt.length; ++i) {
      contigs[i] = new ContigString(contigNt[i]);
    }

    final long[] p1 = {+1, -2};
    final long[] p2 = {+2, -3, +4};
    final long[] p3 = {-3, +5, +6};
    final long[] p4 = {+3, +7};
    final long[][] paths = {p1, p2, p3, p4};

    final GraphImplementation graph = new GraphImplementation(0, Collections.emptyMap(), Collections.emptyMap());
    buildGraph(graph, contigs, paths);

    //graph.dump(System.err);
    checkGraph(graph, contigs, paths, false);
    checkGraph((GraphImplementation) graph.compact(), contigs, paths, false);

    //various bad stuff
    try {
      graph.contig(0);
      fail();
    } catch (final IllegalArgumentException e) {
      //expected
    }
    try {
      graph.pathContig(0, 0);
      fail();
    } catch (final IllegalArgumentException e) {
      //expected
    }
    try {
      graph.pathContig(2, -1);
      fail();
    } catch (final IllegalArgumentException e) {
      //expected
    }
    try {
      graph.pathContig(-2, -1);
      fail();
    } catch (final IllegalArgumentException e) {
      //expected
    }
    try {
      graph.pathContig(2, 3);
      fail();
    } catch (final IllegalArgumentException e) {
      assertEquals("2:3", e.getMessage());
    }
    try {
      graph.pathContig(-2, 3);
      fail();
    } catch (final IllegalArgumentException e) {
      assertEquals("-2:3", e.getMessage());
    }
    try {
      graph.nt(0, 0);
      fail();
    } catch (final IllegalArgumentException e) {
      //expected
    }
    try {
      graph.nt(2, -1);
      fail();
    } catch (final IllegalArgumentException e) {
      assertEquals("2:-1", e.getMessage());
    }
    try {
      graph.nt(-2, -1);
      fail();
    } catch (final IllegalArgumentException e) {
      assertEquals("-2:-1", e.getMessage());
    }
    try {
      graph.nt(2, 4);
      fail();
    } catch (final IllegalArgumentException e) {
      assertEquals("2:4", e.getMessage());
    }
    try {
      graph.nt(-2, 4);
      fail();
    } catch (final IllegalArgumentException e) {
      assertEquals("-2:4", e.getMessage());
    }
    try {
      graph.absPath(0);
      fail();
    } catch (final IllegalArgumentException e) {
      //expected
    }
    try {
      graph.absPath(5);
      fail();
    } catch (final IllegalArgumentException e) {
      //expected
    }
    try {
      graph.absContig(0);
      fail();
    } catch (final IllegalArgumentException e) {
      //expected
    }
    try {
      graph.absContig(9);
      fail();
    } catch (final IllegalArgumentException e) {
      //expected
    }
  }

  //normal case for adds then some deletions and compact at end
  public void testDelete() {
    final String[] contigNt = {"ACGT", "CGTA", "GTAC", "TACG", "AA", "CC", "GG", "TT"};
    final Contig[] contigs = new Contig[contigNt.length];
    for (int i = 0; i < contigNt.length; ++i) {
      contigs[i] = new ContigString(contigNt[i]);
    }

    final long[] p1 = {+1, -2};
    final long[] p2 = {+2, -3, +4};
    final long[] p3 = {-3, +5, +6};
    final long[] p4 = {+3, +7};
    final long[][] paths = {p1, p2, p3, p4};

    final GraphImplementation graph = new GraphImplementation(0, Collections.emptyMap(), Collections.emptyMap());
    buildGraph(graph, contigs, paths);

    //graph.dump(System.err);
    for (int i = 1; i <= contigs.length; ++i) {
      assertFalse(graph.contigDeleted(i));
    }
    for (int i = 1; i <= paths.length; ++i) {
      assertFalse(graph.pathDeleted(i));
    }
    checkGraph(graph, contigs, paths, false);
    //do some deleting
    graph.deleteContig(8);
    graph.deleteContig(2);
    graph.deletePath(4);

    assertTrue(graph.pathDeleted(1));
    assertTrue(graph.pathDeleted(2));
    assertFalse(graph.pathDeleted(3));
    assertTrue(graph.pathDeleted(4));

    assertFalse(graph.contigDeleted(1));
    assertTrue(graph.contigDeleted(2));
    assertFalse(graph.contigDeleted(3));
    assertFalse(graph.contigDeleted(4));
    assertFalse(graph.contigDeleted(5));
    assertFalse(graph.contigDeleted(6));
    assertFalse(graph.contigDeleted(7));
    assertTrue(graph.contigDeleted(8));

    //graph.dump(System.err);

    checkGraph(graph, contigs, paths, true);

    final Contig[] contigsDel = new Contig[contigs.length];
    contigsDel[0] = contigs[0];
    contigsDel[1] = null;
    contigsDel[2] = contigs[2];
    contigsDel[3] = contigs[3];
    contigsDel[4] = contigs[4];
    contigsDel[5] = contigs[5];
    contigsDel[6] = contigs[6];
    contigsDel[7] = null;

    final long[][] pathsDel = new long[paths.length][];
    pathsDel[0] = null;
    pathsDel[1] = null;
    pathsDel[2] = paths[2];
    pathsDel[3] = null;
    checkGraph(graph, contigsDel, pathsDel, false);


    final Contig[] contigsCompact = new Contig[contigs.length - 2];
    contigsCompact[0] = contigs[0];
    contigsCompact[1] = contigs[2];
    contigsCompact[2] = contigs[3];
    contigsCompact[3] = contigs[4];
    contigsCompact[4] = contigs[5];
    contigsCompact[5] = contigs[6];

    final long[][] pathsCompact = new long[1][];
    pathsCompact[0] = new long[] {-2, +4, +5};
    checkGraph((GraphImplementation) graph.compact(), contigsCompact, pathsCompact, false);
  }

  /**
   * @return a handy graph for testing.
   */
  public static GraphImplementation graph() {
    final String[] contigNt = {"ACGT", "CGTA", "GTAC", "TACG", "AA", "CC", "GG", "TT"};
    final Contig[] contigs = new Contig[contigNt.length];
    for (int i = 0; i < contigNt.length; ++i) {
      contigs[i] = new ContigString(contigNt[i]);
    }

    final long[] p1 = {+1, -2};
    final long[] p2 = {+2, -3, +4};
    final long[] p3 = {-3, +5, +6};
    final long[] p4 = {+3, +7};
    final long[][] paths = {p1, p2, p3, p4};

    final Map<String, String> contigAttributes = new HashMap<>();
    contigAttributes.put("foo", "foo comment");
    contigAttributes.put("bar", "bar comment");

    final Map<String, String> pathAttributes = new HashMap<>();
    pathAttributes.put("boo", "comment boo");
    pathAttributes.put("oob", "comment oob");

    final GraphImplementation graph = new GraphImplementation(0, contigAttributes, pathAttributes);
    buildGraph(graph, contigs, paths);
    return graph;
  }

  static GraphImplementation buildGraph(final GraphImplementation graph, final Contig[] contigs, final long[][] paths) {
    graph.globalIntegrity();
    assertEquals(0, graph.numberContigs());
    assertEquals(0, graph.numberPaths());

    for (int i = 0; i < contigs.length; ++i) {
      assertEquals(i + 1, graph.addContig(contigs[i]));
    }

    graph.globalIntegrity();
    assertEquals(contigs.length, graph.numberContigs());
    assertEquals(0, graph.numberPaths());

    for (int i = 0; i < paths.length; ++i) {
      assertEquals(i + 1, graph.addPath(new PathArray(paths[i])));
    }
    return graph;
  }

  private void checkGraph(final GraphImplementation graph, Contig[] contigs, final long[][] paths, final boolean deleted) {
    graph.globalIntegrity();
    if (!deleted) {
      assertEquals(contigs.length, graph.numberContigs());
      assertEquals(paths.length, graph.numberPaths());
    }


    for (int i = 0; i < contigs.length; ++i) {
      checkNt(graph, i + 1, contigs[i], deleted);
    }

    for (int i = 0; i < paths.length; ++i) {
      checkPath(graph, i + 1, paths[i], deleted);
    }

    for (int i = 0; i < contigs.length; ++i) {
      int count = 0;
      for (final long[] path : paths) {
        if (path == null) {
          continue;
        }
        for (final long c : path) {
          if (c == (i + 1) || -c == (i + 1)) {
            ++count;
          }
        }
      }
      if (count == 0) {
        continue;
      }
      final long[] index = new long[count];
      int count1 = 0;
      for (int j = 0; j < paths.length; ++j) {
        final long[] path = paths[j];
        if (path == null) {
          continue;
        }
        for (final long c : path) {
          if (c == (i + 1)) {
            index[count1] = j + 1;
            ++count1;
          } else if (-c == (i + 1)) {
            index[count1] = -(j + 1);
            ++count1;
          }
        }
      }
      checkPaths(graph, i + 1, index, deleted);
    }
  }

  private static long[] neg(long[] paths) {
    final long[] neg = new long[paths.length];
    for (int i = 0; i < paths.length; ++i) {
      final int j = paths.length - i - 1;
      neg[j] = -paths[i];
    }
    return neg;
  }

  private void checkPaths(Graph graph, long contigId, long[] paths, boolean deleted) {
    final long[] sorted = Arrays.copyOf(paths, paths.length);
    Arrays.sort(sorted);
    checkPathsX(graph, contigId, sorted, deleted);
    checkPathsX(graph, -contigId, neg(sorted), deleted);
  }

  private void checkPathsX(Graph graph, long contigId, long[] paths, boolean deleted) {
    //System.err.println("deleted=" + deleted);
    final PathsIterator it = graph.paths(contigId, deleted);
    checkPathsX(graph, contigId, paths, deleted, it);
    if (!deleted) {
      final PathsIterator itd = graph.paths(contigId);
      checkPathsX(graph, contigId, paths, deleted, itd);
    }
    final PathsIterator its = graph.iterator();
    its.set(contigId, deleted);
    checkPathsX(graph, contigId, paths, deleted, its);
  }

  private void checkPathsX(Graph graph, long contigId, long[] paths, boolean deleted, final PathsIterator it) {
    // the order is not guaranteed
    int count = 0;
    final long[] actual = new long[paths.length];

    while (true) {
      final long path = it.nextPathId();
      //System.err.println("contigId=" + contigId + " path=" + path + " length=" + paths.length + " count=" + count); // + " it=" + it.toString());
      if (path == 0) {
        break;
      }
      final int contigIndex = it.contigIndex();
      assertEquals("contigIndex=" + contigIndex, contigId, graph.pathContig(path, contigIndex));
      actual[count] = path;
      ++count;
    }
    Arrays.sort(actual);
    assertEquals(paths.length, count);
    //System.err.println(Arrays.toString(paths));
    //System.err.println(Arrays.toString(actual));
    assertTrue(Arrays.equals(paths, actual));
  }

  private void checkPath(GraphImplementation graph, long pathId, long[] path, boolean deleted) {
    if (path == null) {
      assertTrue(graph.pathDeleted(pathId));
      return;
    }
    final int length = path.length;
    assertEquals(length, graph.pathLength(pathId));
    assertEquals(length, graph.pathLength(-pathId));
    final Path pathp = graph.path(pathId);
    assertEquals(length, pathp.length());
    final Path pathm = graph.path(-pathId);
    assertEquals(length, pathm.length());
    //If any contig has been deleted then the path should be deleted
    boolean contigDeleted = false;
    long max = 0;
    for (int i = 0; i < length; ++i) {
      max = Math.max(max, graph.absContig(path[i]));
      assertEquals(path[i], graph.pathContig(pathId, i));
      assertEquals(path[i], pathp.contig(i));
      assertEquals(i, pathp.index(path[i]));

      final int r = length - i - 1;
      assertEquals(-path[r], graph.pathContig(-pathId, i));
      assertEquals(-path[r], pathm.contig(i));
      assertEquals(i, pathm.index(-path[r]));
      if (graph.contigDeleted(path[i])) {
        contigDeleted = true;
      }
    }
    assertEquals(-1, pathp.index(max + 1));
    assertEquals(-1, pathm.index(max + 1));
    assertEquals(-1, pathp.index(-(max + 1)));
    assertEquals(-1, pathm.index(-(max + 1)));
    if (contigDeleted) {
      assertTrue(graph.pathDeleted(pathId));
    }
  }

  private void checkNt(Graph graph, long contigId, Contig contig, boolean deleted) {
    if (contig == null) {
      assertTrue(graph.contigDeleted(contigId));
      return;
    }
    final int length = contig.length();
    assertEquals(length, graph.contigLength(contigId));
    assertEquals(length, graph.contigLength(-contigId));
    final Contig contigp = graph.contig(contigId);
    Exam.integrity(contigp);
    assertEquals(length, contigp.length());
    final Contig contigm = graph.contig(-contigId);
    Exam.integrity(contigm);
    assertEquals(length, contigm.length());
    for (int i = 0; i < length; ++i) {
      assertEquals("contigId=" + contigId + " i=" + i, contig.nt(i), graph.nt(contigId, i));
      assertEquals(contig.nt(i), contigp.nt(i));

      final byte nt = contig.nt(length - i - 1);
      final byte ntc = DNARange.complement(nt);
      assertEquals(contigId + ":" + i, ntc, graph.nt(-contigId, i));
      assertEquals(contigId + ":" + i, ntc, contigm.nt(i));
    }
    try {
      contigp.nt(contigp.length());
      fail();
    } catch (RuntimeException e) {
      // ok
    }
    try {
      contigp.nt(-1);
      fail();
    } catch (RuntimeException e) {
      // ok
    }
  }
}
