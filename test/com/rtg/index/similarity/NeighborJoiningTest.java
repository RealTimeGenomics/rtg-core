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
package com.rtg.index.similarity;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import com.rtg.mode.SequenceType;
import com.rtg.reader.MockArraySequencesReader;
import com.rtg.reader.SequencesReader;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.MemoryPrintStream;

import junit.framework.TestCase;

/**
 * Tests the corresponding class.
 *
 */
public class NeighborJoiningTest extends TestCase {

  private static final String LS = System.lineSeparator();

  private static final HashMap<String, Integer> VALID = new HashMap<>();
  static {
    VALID.put("X   D" + LS
        + "X  4" + LS
        + "X   C" + LS
        + "X 3" + LS
        + "X  B" + LS
        + "X2" + LS
        + "X A" + LS, 0);
    VALID.put("X   D" + LS
        + "X  4" + LS
        + "X   C" + LS
        + "X 3" + LS
        + "X  A" + LS
        + "X2" + LS
        + "X B" + LS, 1);
    VALID.put("X  D" + LS
        + "X 3" + LS
        + "X  C" + LS
        + "X2" + LS
        + "X  B" + LS
        + "X 4" + LS
        + "X  A" + LS, 2);
    VALID.put("X   B" + LS
        + "X  4" + LS
        + "X   A" + LS
        + "X 3" + LS
        + "X  C" + LS
        + "X2" + LS
        + "X D" + LS, 3);
    VALID.put("X  B" + LS
        + "X 3" + LS
        + "X  A" + LS
        + "X2" + LS
        + "X  D" + LS
        + "X 4" + LS
        + "X  C" + LS, 4);
    VALID.put("X   B" + LS
        + "X  4" + LS
        + "X   A" + LS
        + "X 3" + LS
        + "X  D" + LS
        + "X2" + LS
        + "X C" + LS, 5);
  }

  static ArrayList<ArrayList<Double>> buildList(final double[][] d) {
    final ArrayList<ArrayList<Double>> res = new ArrayList<>();
    for (final double[] dd : d) {
      final ArrayList<Double> row = new ArrayList<>();
      for (double aDd : dd) {
        row.add(aDd);
      }
      res.add(row);
    }
    return res;
  }

  @Override
  public void setUp() {
    Diagnostic.setLogStream();
  }

  @Override
  public void tearDown() {
    Diagnostic.setLogStream();
  }

  public void testWikiExample() throws IOException {
    final MemoryPrintStream log = new MemoryPrintStream();
    Diagnostic.setLogStream(log.printStream());
    final double[][] d = {
        {},
        {7},
        {11, 6},
        {14, 9, 7},
    };

    final ArrayList<String> s = new ArrayList<>();
    s.add("A");
    s.add("B");
    s.add("C");
    s.add("D");
    final int[] counts = new int[VALID.size()];
    //the result is non-deterministic - do enough to get a statistical result.
    for (int k = 0; k < 100; k++) {
      final NeighborJoining nj = new NeighborJoining(k * 1237L);
      final BinaryTree tree = nj.neighborJoin(s, buildList(d));
      assertNotNull(tree);
      final ByteArrayOutputStream bos = new ByteArrayOutputStream();
      try {
        try (PrintStream ps = new PrintStream(bos)) {
          final double dist = tree.dump("X", ps);
          assertTrue("" + dist, dist > 16.9);
          assertTrue("" + dist, dist < 17.1);
        }
      } finally {
        bos.close();
      }
      final String res = bos.toString();
      final Integer v = VALID.get(res);
      if (v == null) {
        fail(res);
      } else {
        counts[v]++;
      }
    }
    for (final int count : counts) {
      assertTrue(count > 2);
    }
    assertTrue(log.toString().contains("NeighborJoining: B + A -> 4"));
  }

  public void testWikiExample2() throws IOException {
    final double[][] d = {
        {},
        {7},
        {11, 6},
        {14, 9, 7},
    };

    final ArrayList<String> s = new ArrayList<>();
    s.add("A");
    s.add("B");
    s.add("C");
    s.add("D");
    final NeighborJoining nj = new NeighborJoining();
    final BinaryTree tree = nj.neighborJoin(s, buildList(d));
    assertNotNull(tree);
    final Appendable out = new StringWriter();
    tree.modifiedNewHampshire(out);
    final String nh = out.toString();
    final ByteArrayOutputStream bos = new ByteArrayOutputStream();
    try {
      try (final PrintStream ps = new PrintStream(bos)) {
        final double dist = tree.dump("X", ps);
        assertTrue("" + dist, dist > 16.9);
        assertTrue("" + dist, dist < 17.1);
      }
    } finally {
      bos.close();
    }
    final String dmp = bos.toString();
//  System.err.println(nh);
//  System.err.println(dmp);

    final String expnh = ""
      + "(" + LS
      + " (" + LS
      + "  (" + LS
      + "   B:1.0000," + LS
      + "   A:6.0000" + LS
      + "  ):3.0000," + LS
      + "  D:5.0000" + LS
      + " ):1.0000," + LS
      + " C:1.0000" + LS
      + ")" + LS
      ;
    assertEquals(expnh, nh);

    final Appendable xmlOut = new StringWriter();
    tree.phyloXml(xmlOut);
    final String xml = xmlOut.toString();
    assertEquals("<?xml version=\"1.0\" encoding=\"UTF-8\"?>" + LS
                 + "<phyloxml xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://www.phyloxml.org http://www.phyloxml.org/1.10/phyloxml.xsd\" xmlns=\"http://www.phyloxml.org\">" + LS
                 + "  <phylogeny rooted=\"false\">" + LS
                 + "    <description>Tree produced by the RTG phylogeny module</description>" + LS
                 + "    <clade>" + LS
                 + "      <clade branch_length=\"1.0000\">" + LS
                 + "        <clade branch_length=\"3.0000\">" + LS
                 + "          <clade branch_length=\"1.0000\">" + LS
                 + "            <name>B</name>" + LS
                 + "          </clade>" + LS
                 + "          <clade branch_length=\"6.0000\">" + LS
                 + "            <name>A</name>" + LS
                 + "          </clade>" + LS
                 + "        </clade>" + LS
                 + "        <clade branch_length=\"5.0000\">" + LS
                 + "          <name>D</name>" + LS
                 + "        </clade>" + LS
                 + "      </clade>" + LS
                 + "      <clade branch_length=\"1.0000\">" + LS
                 + "        <name>C</name>" + LS
                 + "      </clade>" + LS
                 + "    </clade>" + LS
                 + "  </phylogeny>" + LS
                 + "</phyloxml>" + LS, xml);

    final String expdmp = ""
      + "X   B" + LS
      + "X  4" + LS
      + "X   A" + LS
      + "X 3" + LS
      + "X  D" + LS
      + "X2" + LS
      + "X C" + LS
      ;
    assertEquals(expdmp, dmp);

    final String expstr = ""
      + "node 0 \"B\"" + LS
      + "node 1 \"4\"" + LS
      + "node 2 \"A\"" + LS
      + "node 3 \"3\"" + LS
      + "node 4 \"D\"" + LS
      + "node 5 \"2\"" + LS
      + "node 6 \"C\"" + LS
      + "edge 5-3" + LS
      + "edge 3-1" + LS
      + "edge 1-0" + LS
      + "edge 1-2" + LS
      + "edge 3-4" + LS
      + "edge 5-6" + LS
      ;
    assertEquals(expstr, tree.toString());
  }

  public void testMakeMatrix() {
    final SimilarityMatrix ma = new SimilarityMatrix(4);
    final int[][] n = {
        {1},
        {0, 4},
        {2, 2, 9},
        {2, 4, 2, 16}
    };
    for (int i = 0; i < n.length; i++) {
      for (int j = 0; j < n[i].length; j++) {
        ma.set(i, j, n[i][j]);
      }
    }
    final StringWriter sb = new StringWriter();
    final ArrayList<ArrayList<Double>> aa = NeighborJoining.makeArray(ma);
    for (ArrayList<Double> a : aa) {
      for (Double d : a) {
        final String strVal = String.format("  %1$02.4f", d);
        sb.append(strVal);
      }
      sb.append(LS);
    }

    //System.err.println(sb.toString());
    final String exp = ""
      + "" + LS
      + "  1.0000" + LS
      + "  0.6000  0.7500" + LS
      + "  0.6667  0.6667  0.8571" + LS
      ;
    assertEquals(exp, sb.toString());
  }

  public void testNeighborJoin1() throws IOException {
    final SequencesReader reader = new MockArraySequencesReader(SequenceType.DNA, 4);
    final ArrayList<String> na = makeNames(reader);
    assertEquals("[seq0, seq1, seq2, seq3]", Arrays.toString(na.toArray()));

    final SimilarityMatrix ma = new SimilarityMatrix(4);
    final int[][] n = {
        {1},
        {0, 4},
        {2, 2, 9},
        {2, 4, 2, 16}
    };
    for (int i = 0; i < n.length; i++) {
      for (int j = 0; j < n[i].length; j++) {
        ma.set(i, j, n[i][j]);
      }
    }

    final NeighborJoining nj = new NeighborJoining();

    final BinaryTree tree = nj.neighborJoin(na, ma);
    final Appendable out = new StringWriter();
    tree.modifiedNewHampshire(out);
    final String str = out.toString();
    //System.err.println(str);
    final String exp = ""
      + "(" + LS
      + " (" + LS
      + "  seq3:0.2768," + LS
      + "  seq1:0.3899" + LS
      + " ):0.0926," + LS
      + " (" + LS
      + "  seq2:0.2851," + LS
      + "  seq0:0.3149" + LS
      + " ):0.0926" + LS
      + ")" + LS
      ;
    assertEquals(exp, str);
  }

  /** A bug seen in the wild. */
  public void testNeighborJoin2() throws IOException {
    final SequencesReader reader = new MockArraySequencesReader(SequenceType.DNA, 3);
    final ArrayList<String> na = makeNames(reader);
    assertEquals("[seq0, seq1, seq2]", Arrays.toString(na.toArray()));

    final SimilarityMatrix ma = new SimilarityMatrix(3);
    final int[][] n = {
        {0},
        {0, 2},
        {0, 0, 0},
    };
    for (int i = 0; i < n.length; i++) {
      for (int j = 0; j < n[i].length; j++) {
        ma.set(i, j, n[i][j]);
      }
    }

    final NeighborJoining nj = new NeighborJoining();

    final BinaryTree tree = nj.neighborJoin(na, ma);
    final Appendable out = new StringWriter();
    tree.modifiedNewHampshire(out);
    final String str = out.toString();
    //System.err.println(str);
    final String exp = ""
      + "(" + LS
      + " (" + LS
      + "  seq1:0.5000," + LS
      + "  seq0:0.5000" + LS
      + " ):0.2500," + LS
      + " seq2:0.2500" + LS
      + ")" + LS
      ;
    assertEquals(exp, str);
  }

  /**
   * Get all sequence names in the reader.
   * @param reader sequences from where to get names.
   * @return list of all sequence names in the reader.
   */
  static ArrayList<String> makeNames(final SequencesReader reader) throws IOException {
    final ArrayList<String> nodeNames;
    nodeNames = new ArrayList<>();
    for (int i = 0; i < reader.numberSequences(); i++) {
      reader.seek(i);
      nodeNames.add(reader.currentName());
    }
    return nodeNames;
  }

}
