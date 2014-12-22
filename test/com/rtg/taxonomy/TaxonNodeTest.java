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

import java.util.HashMap;
import java.util.List;

import junit.framework.TestCase;

/**
 */
public class TaxonNodeTest extends TestCase {

  public void testSetters() {
    final TaxonNode tn = new TaxonNode(123);
    assertEquals(123, tn.getId());
    assertNull(tn.getName());
    assertNull(tn.getRank());
    assertNull(tn.getParent());
    assertEquals(-1, tn.getParentId());

    tn.setName("name");
    tn.setRank("rank");

    assertEquals(123, tn.getId());
    assertEquals("name", tn.getName());
    assertEquals("rank", tn.getRank());
    assertNull(tn.getParent());
    assertEquals(-1, tn.getParentId());

    final TaxonNode tn2 = new TaxonNode(123, "name", "rank");
    assertTrue(tn.equals(tn2));
    assertEquals(0, tn.compareTo(tn2));

    tn2.setName("name2");
    tn2.setRank("rank\ta bc\n");

    assertTrue(tn.equals(tn2));
    assertEquals(0, tn.compareTo(tn2));

    assertEquals(123, tn2.getId());
    assertEquals("name2", tn2.getName());
    assertEquals("rank a bc ", tn2.getRank());
    assertNull(tn2.getParent());
    assertEquals(-1, tn2.getParentId());
  }

  public void testSimpleTree() {
    final HashMap<Integer, TaxonNode> nodes = new HashMap<>();
    final TaxonNode root = new TaxonNode(1, "root", "root");
    nodes.put(1, root);
    for (int i = 2; i <= 10; i++) {
      TaxonNode node = new TaxonNode(i, "node" + i, "rank" + i);
      nodes.put(i, node);
      nodes.get((i - 2) / 3 + 1).addChild(node);
    }

    List<TaxonNode> traverse = root.depthFirstTraversal();
    assertEquals(10, traverse.size());

    int i = 0;
    int[] ids = {1, 2, 5, 6, 7, 3, 8, 9, 10, 4};
    for (TaxonNode tn : traverse) {
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

    final TaxonNode node = nodes.get(10);
    node.detach();
    traverse = root.depthFirstTraversal();
    assertEquals(9, traverse.size());
  }

  public void testCommon() {
    final TaxonNode tn = new TaxonNode(123, "name", "rank");
    assertEquals("123\t-1\trank\tname", tn.toString());
    assertEquals(123, tn.hashCode());
    assertFalse(tn.equals(null));
    assertTrue(tn.equals(tn));
    assertEquals(0, tn.compareTo(tn));

    final TaxonNode tn2 = new TaxonNode(234, "name2", "rank2");
    assertEquals(111, tn2.compareTo(tn));
    assertEquals(-111, tn.compareTo(tn2));
    assertFalse(tn.equals(tn2));
  }
}
