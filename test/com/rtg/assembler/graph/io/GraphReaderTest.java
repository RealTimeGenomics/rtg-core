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
package com.rtg.assembler.graph.io;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.UUID;

import com.rtg.assembler.GraphMapCliTest;
import com.rtg.assembler.graph.Contig;
import com.rtg.assembler.graph.Graph;
import com.rtg.assembler.graph.GraphFactory;
import com.rtg.assembler.graph.MutableGraph;
import com.rtg.assembler.graph.Path;
import com.rtg.util.Pair;
import com.rtg.util.store.StoreDirProxy;
import com.rtg.util.store.StoreDirResource;
import com.rtg.util.store.StoreDirString;
import com.rtg.util.store.StoreDirectory;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 */
public class GraphReaderTest extends TestCase {

  private static final String RESOURCE = "com/rtg/assembler/graph/io/resources";

  public void testTopLevelBogusInput() {
    try {
      GraphReader.read(new StoreDirProxy(new File("what-a-stupid-name-for-a-file")));
      fail();
    } catch (final IOException e) {
      final String message = e.getMessage();
      assertEquals("Not a directory what-a-stupid-name-for-a-file", message);
    }
  }

  public void testMissingManifest() throws IOException {
    final File root = FileHelper.createTempDirectory();
    try {
      FileHelper.resourceToFile(RESOURCE + "/" + "header1", new File(root, GraphReader.HEADER_LEGACY));
      final StoreDirectory rootDir = new StoreDirProxy(root);
      try {
        GraphReader.read(rootDir);
        fail();
      } catch (final IOException e) {
        assertTrue(e.getMessage().contains("Missing manifest.txt"));
      }
    } finally {
      FileHelper.deleteAll(root);
    }
  }

  public void testMissingHeader() throws IOException {
    final File root = FileHelper.createTempDirectory();
    try {
      FileHelper.resourceToFile(RESOURCE + "/" + "manifest", new File(root, GraphReader.MANIFEST));
      final StoreDirectory rootDir = new StoreDirProxy(root);
      try {
        GraphReader.read(rootDir);
        fail();
      } catch (final IOException e) {
        assertTrue(e.getMessage().contains("Missing header.tsv"));
      }
    } finally {
      FileHelper.deleteAll(root);
    }
  }

  public void testKnownGraph() throws IOException {
    final StoreDirectory rootDir = new StoreDirResource("com/rtg/assembler/planning/resources");
    final Graph g = GraphReader.read(rootDir);
    checkGraph(g);
  }
  public void testKnownGraphPlusAttributes() throws IOException {
    final Map<String, String> contigAttrs = new HashMap<>();
    contigAttrs.put("hairiness", "how much hair does this contig have");
    final Map<String, String> pathAttrs = new HashMap<>();
    pathAttrs.put("BMI", "rough proxy for how fat the path is");
    final StoreDirectory rootDir = new StoreDirResource("com/rtg/assembler/planning/resources");
    final Graph g = GraphReader.read(GraphFactory.KMER, rootDir, contigAttrs, pathAttrs);
    assertTrue(g.contigAttributes().containsKey("hairiness"));
    assertTrue(g.pathAttributes().containsKey("BMI"));
  }

  public void testKnownGraphDefault() throws IOException {
    final StoreDirectory rootDir = new StoreDirResource("com/rtg/assembler/planning/resources");
    final Graph g = GraphReader.read(GraphFactory.DEFAULT, rootDir);
    checkGraph(g);
  }

  public void testBadManifestIndex() {
    checkBadFile("manifest.txt invalid line:  header.txt", "762b36734794180251a4a46aefde8a68", "");
  }

  public void testBadPathContigLine() {
    checkBadFile("Invalid contig identifier: 242 in line:contigs\t+242\t-4", "\t-1", "42");
  }

  public void testBadPathTooFewContigs() {
    checkBadFile("Insufficient (<2) contig identifiers on path: contigs\t-2", "\\+1\t", "");
  }

  public void testBadContigHeader() {
    checkBadFile("Invalid contig header line:>\tcov=42", ">1\t", ">\t");
  }

  public void testBadInputGuidHeader() {
    checkBadFile("header.txt: bad attribute line: foo\t3F2504E0-4F89-11D3-9A0C-0305E82C3302", "inputguid", "foo");
  }

  //too many fields
  public void testBadHeaderAttribute() {
    checkBadFile("header.txt: bad attribute line: contig\tcov\tfoo\tCoverage of the contig", "contig\tcov", "contig\tcov\tfoo");
  }

  public void testBadHeaderVersion() {
    checkBadFile("header.txt: expected version line, saw: version", "version\t0.0 head", "version");
  }

  public void testBadHeaderDate() {
    checkBadFile("header.txt: malformed date line, saw: foo\t2012/05/07 14:53:05\t", "date", "foo");
  }

  //nothing after version
  public void testBadHeader() {
    final Set<Pair<String, String>> rules = new HashSet<>();
    rules.add(new Pair<>("date", "#date"));
    rules.add(new Pair<>("command", "#command"));
    rules.add(new Pair<>("guid", "#guid"));
    rules.add(new Pair<>("inputguid", "#inputguid"));
    rules.add(new Pair<>("contig", "#contig"));
    rules.add(new Pair<>("path", "#path"));
    checkBadFile("header.txt: missing or incorrect date line", rules);
  }

  //too few fields fields
  public void testBadHeaderTooFew0() {
    checkBadFile("header.txt: expected guid line, saw: guid", "guid\t3F2504E0-4F89-11D3-9A0C-0305E82C3301", "guid");
  }

  //too few fields fields
  public void testBadHeaderTooFew1() {
    checkBadFile("header.txt: bad attribute line: contig", "contig\tfoo", "contig");
  }

  public void testBadContigToken() {
    checkBadFile("Expected line starting with 'contigs', saw: foo\t+2\t-1\t-4", "contigs\t\\+2", "foo\t\\+2");
  }

  public void testBadContigTokenLast() {
    checkBadFile("Missing contigs after: path\t2", "contigs\t\\+1", "#contigs\\t\\\\+1");
  }

  public void testBadPathTokenFirst() {
    checkBadFile("Expected line starting with 'path', saw: foo\t1\tmin=3", "path\t1", "foo\t1");
  }

  public void testBadManifest() {
    checkBadFile("header.txt: md5sum error, expected 73dd4e9a93da42ca682b0868edf3f42d got 762b36734794180251a4a46aefde8a68", "762b36734794180251a4a46aefde8a68", "73dd4e9a93da42ca682b0868edf3f42d");
  }

  public void testBadNoDate() {
    checkBadFile("header.txt: expected date line, saw: foobar2012/05/07 14:53:05\t", "date\t", "foobar");
  }

  private void checkBadFile(final String msg, final String from, final String to) {
    final Set<Pair<String, String>> rules = new HashSet<>();
    rules.add(new Pair<>(from, to));
    checkBadFile(msg, rules);
  }

  private void checkBadFile(final String msg, final Set<Pair<String, String>> rules) {
    final StoreDirectory rootDir = new StoreDirResource("com/rtg/assembler/planning/resources");
    final StoreDirectory mutant = new StoreDirMutator(rootDir, rules);
    try {
      GraphReader.read(mutant);
      fail();
    } catch (final Exception e) {
      //e.printStackTrace();
      assertEquals(msg, e.getMessage());
    }
  }

  static void checkGraph(final Graph g) {
    assertEquals(5, g.numberContigs());
    assertEquals(2, g.numberPaths());
    final Map<String, String> ca = g.contigAttributes();
    assertEquals(4, ca.size());
    assertEquals("Coverage of the contig", ca.get("cov"));
    assertEquals("Tip value", ca.get("tip"));
    assertEquals("", ca.get("foo"));

    final Map<String, String> pa = g.pathAttributes();
    assertEquals(2, pa.size());
    assertTrue(pa.containsKey("min"));
    // Following depends on sequential nature of load an id allocation in graph
    assertEquals(10, g.contigLength(1));
    final Contig c1 = g.contig(1);
    assertEquals(1, c1.nt(0));
    assertEquals(2, c1.nt(1));
    assertEquals(3, c1.nt(2));
    assertEquals(4, c1.nt(3));
    assertEquals(0, c1.nt(4));
    assertEquals(0, c1.nt(5));
    assertEquals(0, c1.nt(6));
    assertEquals(0, c1.nt(7));
    assertEquals(0, c1.nt(8));
    assertEquals(0, c1.nt(9));
    assertEquals(1, g.nt(1, 0));
    assertEquals(4, g.nt(1, 3));
    assertEquals(1, g.nt(3, 0));
    assertEquals(4, g.nt(3, 3));
    assertEquals("AAAAAAAAAAAAAAA".length(), g.contigLength(2));
    assertEquals("AAAAAAAAAAAAAAA".length(), g.contigLength(4));
    assertEquals("42", g.contigAttribute(3, "cov"));
    assertEquals("1.9", g.contigAttribute(3, "tip"));
    assertEquals("12.34", g.contigAttribute(2, "tip"));
    assertEquals("1", g.contigAttribute(4, "cov"));
    assertNull(g.pathAttribute(2, "min"));
    assertEquals("3", g.pathAttribute(1, "min"));
    final Path p1 = g.path(1);
    assertEquals(3, p1.length());
    assertEquals(2, p1.contig(0));
    assertEquals(-1, p1.contig(1));
    assertEquals(-4, p1.contig(2));
    final Path p2 = g.path(2);
    assertEquals(2, p2.length());
    assertEquals(+1, p2.contig(0));
    assertEquals(-2, p2.contig(1));
  }

  public void testGuid() throws IOException {

    final MutableGraph builtGraph = GraphMapCliTest.makeGraph(2
        , new String[]{"ACGT", "GGGG", "TTAA"}
        , new long[][]{{1, 2}, {3, 2}});
    StoreDirString storeDir = new StoreDirString();
    GraphWriter.write(builtGraph, storeDir, "monkey", Collections.<UUID>emptySet());
    GraphReader.getUUID(storeDir);
    storeDir = new StoreDirString();
    try (OutputStream out = storeDir.child("header.txt").outputStream()) {
      out.write("nananananananana".getBytes());
    }
    try {
      GraphReader.getUUID(storeDir);
      fail();
    } catch (IOException e) {
      assertEquals("UUID for graph is missing", e.getMessage());
    }

  }
}
