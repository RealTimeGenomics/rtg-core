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

package com.rtg.taxonomy;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.StringWriter;
import java.util.ArrayList;

import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;
import com.rtg.util.test.NanoRegression;

import junit.framework.TestCase;

/**
 * Tests for Taxonomy class.
 */
public class TaxonomyTest extends TestCase {

  protected NanoRegression mNano = null;

  @Override
  public void setUp() {
    Diagnostic.setLogStream();
    mNano = new NanoRegression(this.getClass(), false);
  }

  @Override
  public void tearDown() throws Exception {
    // clear the module name so later tests don't report SlimException to the
    // Talkback system
    try {
      mNano.finish();
    } finally {
      mNano = null;
    }
  }

  public void testConstructor() {
    final Taxonomy tax = new Taxonomy();
    assertEquals(0, tax.size());
    assertNull(tax.getRoot());
    assertEquals(0, tax.subset(new ArrayList<Integer>()).size());

    assertEquals("taxonomy.tsv", TaxonomyUtils.TAXONOMY_FILE);
    assertEquals(1, Taxonomy.ROOT_ID);
  }

  public void testAddNode() {
    final Taxonomy tax = new Taxonomy();

    tax.addNode(1, -1, "root", "root");
    for (int i = 2; i <= 10; i++) {
      tax.addNode(i, (i - 2) / 3 + 1, "node" + i, "rank" + i);
    }
    assertEquals(10, tax.size());
    assertNotNull(tax.getRoot());

    assertTrue(tax.isConsistent());

    int i = 0;
    int[] ids = {1, 2, 5, 6, 7, 3, 8, 9, 10, 4};
    for (TaxonNode tn : tax.getRoot().depthFirstTraversal()) {
      assertEquals(ids[i], tn.getId());
      if (tn.getId() != 1) {
        assertEquals((tn.getId() - 2) / 3 + 1, tn.getParentId());
        assertEquals(tn.getParent().getId(), tn.getParentId());
      } else {
        assertEquals(-1, tn.getParentId());
        assertNull(tn.getParent());
      }
      i++;
    }
    assertTrue(tax.isConsistent());
  }

  public void testRemoveSubTree() {
    final Taxonomy tax = new Taxonomy();

    tax.addNode(1, -1, "root", "root");
    for (int i = 2; i <= 10; i++) {
      tax.addNode(i, (i - 2) / 3 + 1, "node" + i, "rank" + i);
    }
    assertEquals(10, tax.size());
    assertNotNull(tax.getRoot());
    assertTrue(tax.isConsistent());

    assertTrue(tax.removeSubTree(3));
    assertEquals(6, tax.size());
    assertFalse(tax.contains(3));
    assertFalse(tax.contains(10));
    assertFalse(tax.contains(9));
    assertFalse(tax.contains(8));
    assertTrue(tax.contains(7));

    assertFalse(tax.removeSubTree(3));
    assertEquals(6, tax.size());
  }

  public void testSubset() {
    final Taxonomy tax = new Taxonomy();

    tax.addNode(1, -1, "root", "root");
    for (int i = 2; i <= 10; i++) {
      tax.addNode(i, (i - 2) / 3 + 1, "node" + i, "rank" + i);
    }

    final ArrayList<Integer> ids = new ArrayList<>();
    ids.add(2);
    final Taxonomy sub = tax.subset(ids);
    assertEquals(2, sub.size());

    ids.add(11);
    try {
      tax.subset(ids);
      fail("Accepted bad tax id.");
    } catch (IllegalArgumentException e) {
      assertEquals("Taxonomy does not contain node with id 11", e.getMessage());
    }
  }

  public void testReadWrite() throws IOException {
    final File dir = FileUtils.createTempDir("taxonomy", "test");
    try {
      final File tree1 = new File(dir, "tree1.tsv");
      FileHelper.resourceToFile("com/rtg/taxonomy/resources/tree1.tsv", tree1);

      final Taxonomy tax1 = new Taxonomy();
      try (final FileInputStream fis = new FileInputStream(tree1)) {
        tax1.read(fis);
      }
      assertEquals(100, tax1.size());
      assertTrue(tax1.isConsistent());
      assertTrue(tax1.contains(1));

      final ArrayList<Integer> ids = new ArrayList<>();
      for (int id = 156657; id <= 156666; id++) {
        ids.add(id);
        assertTrue(tax1.contains(id));
      }
      assertEquals(10, ids.size());
      final Taxonomy tax2 = tax1.subset(ids);
      assertEquals(14, tax2.size());
      assertTrue(tax2.isConsistent());

      final StringWriter sw = new StringWriter();
      tax2.write(sw);
      mNano.check("tree2.tsv", sw.toString(), false);

    } finally {
      FileHelper.deleteAll(dir);
    }
  }

  public void testConsistent() {
    final Taxonomy tax = new Taxonomy();
    assertFalse(tax.isConsistent());
    assertEquals("No root node.", tax.getInconsistencyReason());

    tax.addNode(1, -1, "root", null);
    assertFalse(tax.isConsistent());
    assertEquals("Node 1 does not have a rank.", tax.getInconsistencyReason());

    tax.clear();
    tax.addNode(1, -1, null, "root");
    assertFalse(tax.isConsistent());
    assertEquals("Node 1 does not have a name.", tax.getInconsistencyReason());

    tax.clear();
    tax.addNode(1, -1, "root", "root");
    assertTrue(tax.isConsistent());

    for (int i = 10; i >= 3; i--) {
      int parId = (i - 2) / 3 + 1;
      tax.addNode(i, parId, "node" + i, "rank" + i);
      assertFalse("i:" + i, tax.isConsistent());
      assertTrue("i:" + i, tax.getInconsistencyReason().contains("does not have a name."));
    }

    tax.addNode(2, 1, "node2", "rank2");
    assertTrue(tax.isConsistent());

    assertEquals(10, tax.size());
    assertNotNull(tax.getRoot());

    int i = 0;
    int[] ids = {1, 2, 5, 6, 7, 3, 8, 9, 10, 4};
    for (TaxonNode tn : tax.getRoot().depthFirstTraversal()) {
      assertEquals(ids[i], tn.getId());
      if (tn.getId() != 1) {
        assertEquals((tn.getId() - 2) / 3 + 1, tn.getParentId());
        assertEquals(tn.getParent().getId(), tn.getParentId());
      } else {
        assertEquals(-1, tn.getParentId());
        assertNull(tn.getParent());
      }
      i++;
    }
    assertTrue(tax.isConsistent());

    tax.addMergeNode(5, 20);
    assertTrue(tax.isConsistent());

    final TaxonNode node = tax.get(20);
    assertNotNull(node);
    assertEquals("same", node.getName());
    assertEquals("same", node.getRank());
    assertEquals(20, node.getId());
    assertEquals(5, node.getParentId());

    tax.addNode(30, -1, "bad", "bad");
    assertFalse(tax.isConsistent());
    assertEquals("Node 30 does not link to root.", tax.getInconsistencyReason());

  }

  private FileInputStream createFileStream(File dir, String content) throws IOException {
    final File file = FileHelper.createTempFile(dir);
    FileUtils.stringToFile(content, file);
    return new FileInputStream(file);
  }

  public void testBadFileFormat() throws Exception {
    final String header = "#RTG taxonomy version 1.0" + StringUtils.LS;
    final String badHeader = "#YYY taxonomy version 1.0" + StringUtils.LS;
    final String badHeader2 = "#RTG taxonomy version YYY" + StringUtils.LS;
    final String comment = "#taxID\tparentID\trank\tname" + StringUtils.LS;
    final String lines = "1\t-1\tno rank\troot" + StringUtils.LS
        + "10239\t1\tsuperkingdom\tViruses" + StringUtils.LS;

    final File dir = FileUtils.createTempDir("taxonomy", "test");
    try {
      final Taxonomy tax = new Taxonomy();
      tax.read(createFileStream(dir, header + comment + lines));
      assertEquals(2, tax.size());
      assertTrue(tax.isConsistent());
      assertTrue(tax.contains(1));
      tax.clear();

      try {
        tax.read(createFileStream(dir, badHeader + comment + lines));
        fail("accepted bad header");
      } catch (IOException e) {
        assertEquals("No version information on first line of file.", e.getMessage());
      }
      tax.clear();

      try {
        tax.read(createFileStream(dir, badHeader2 + comment + lines));
        fail("accepted bad header");
      } catch (IOException e) {
        assertEquals("Expecting version 1.0 but saw YYY", e.getMessage());
      }
      tax.clear();

      for (String line : new String[] {
          "12333\t10239\tno rank unclassified phages",
          "12333\t10239\tno rank\tunclassified\tphages",
          "some text that is not formatted"
      }) {
        try {
          tax.read(createFileStream(dir, header + comment + lines + line));
          fail("accepted bad line: " + line);
        } catch (IOException e) {
          assertEquals("Malformed taxonomy file line 5: " + line, e.getMessage());
        }
        tax.clear();
      }

      try {
        tax.read(createFileStream(dir, header + comment + lines + "12333\t10239\tno rank\tunclassified phages" + StringUtils.LS + "12333\t10239\tno rank\tunclassified phages" + StringUtils.LS));
        fail("accepted duplicate id");
      } catch (IOException e) {
        assertEquals("Duplicate taxon id: 12333", e.getMessage());
      }
      tax.clear();

      final String line = "blah\t10239\tno rank\tunclassified phages";
      try {
        tax.read(createFileStream(dir, header + comment + lines + line + StringUtils.LS));
        fail("bad number");
      } catch (IOException e) {
        assertEquals("Malformed taxonomy file line 5: " + line, e.getMessage());
      }
      tax.clear();


    } finally {
      FileHelper.deleteAll(dir);
    }
  }

}
