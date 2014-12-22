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
package com.rtg.util;

import java.util.Map;
import java.util.Set;

import junit.framework.TestCase;

/**
 */
public class MapSetTest extends TestCase {

  /**
   * Test method for {@link MapSet}.
   */
  public final void test1() {
    final MapSet<String, Integer> ms = new MapSet<>();
    assertEquals(0, ms.entrySet().size());
    assertEquals(0, ms.numberOfKeys());
    assertEquals(0, ms.keySet().size());
    assertFalse(ms.contains("a", 1));
    assertEquals(null, ms.get("a"));
    assertEquals(null, ms.get("b"));
    final int hash0 = ms.hashCode();
    assertTrue(hash0 == 31); //regression test

    ms.put("a", 1);
    assertTrue(ms.contains("a", 1));
    assertFalse(ms.contains("b", 2));
    final Set<Integer> sa = ms.get("a");
    assertEquals(1, sa.size());
    assertEquals(1, ms.numberOfKeys());
    assertEquals(1, ms.keySet().size());
    assertTrue(ms.keySet().contains("a"));
    sa.contains(1);
    assertEquals(null, ms.get("b"));
    final int hash1 = ms.hashCode();
    assertTrue(hash0 != hash1); //low probability this will fail

    ms.put("b", 2);
    assertTrue(ms.contains("a", 1));
    assertTrue(ms.contains("b", 2));
    final Set<Integer> sb = ms.get("b");
    assertEquals(1, sb.size());
    assertEquals(2, ms.numberOfKeys());
    assertEquals(2, ms.keySet().size());
    assertTrue(ms.keySet().contains("a"));
    assertTrue(ms.keySet().contains("b"));
    sb.contains(2);
    final int hash2 = ms.hashCode();
    assertTrue(hash0 != hash1 && hash1 != hash2); //low probability this will fail

    ms.put("a", 1);
    assertEquals(2, ms.numberOfKeys());
    ms.put("a", 3);
    assertEquals(2, ms.numberOfKeys());
    assertEquals(2, ms.keySet().size());
    assertTrue(ms.keySet().contains("a"));
    assertTrue(ms.keySet().contains("b"));
    assertTrue(ms.contains("a", 1));
    assertTrue(ms.contains("b", 2));
    assertTrue(ms.contains("a", 3));
    final Set<Integer> sa2 = ms.get("a");
    assertEquals(2, sa2.size());
    sa.contains(1);
    sa.contains(3);

    final Set<Map.Entry<String, Set<Integer>>> me = ms.entrySet();
    assertEquals(2, me.size());
    for (final Map.Entry<String, Set<Integer>> entry : me) {
      final String key = entry.getKey();
      final Set<Integer> values = entry.getValue();
      if (key.equals("a")) {
        assertEquals(2, values.size());
        assertTrue(values.contains(1));
        assertTrue(values.contains(3));
      } else if (key.equals("b")) {
        assertEquals(1, values.size());
        assertTrue(values.contains(2));
      } else {
        fail();
      }
    }

    assertEquals(sa, ms.remove("a"));
    assertFalse(ms.contains("a", 1));
    assertTrue(ms.contains("b", 2));
    final Set<Integer> sb1 = ms.get("b");
    assertEquals(1, sb1.size());
    assertEquals(1, ms.numberOfKeys());
    assertEquals(1, ms.keySet().size());
    assertFalse(ms.keySet().contains("a"));
    assertTrue(ms.keySet().contains("b"));
    assertEquals(null, ms.get("a"));

    assertEquals(null, ms.remove("a"));

  }

  public void testEquals1() {
    final MapSet<String, Integer> ms0 = new MapSet<>();
    final MapSet<String, Integer> msa = new MapSet<>();
    msa.put("a", 1);
    msa.put("b", 2);
    msa.put("a", 1);
    final MapSet<String, Integer> msb = new MapSet<>();
    msb.put("a", 1);
    msb.put("b", 2);
    msb.put("a", 1);
    TestUtils.equalsHashTest(new MapSet[][] {{ms0}, {msa, msb}});
    assertFalse(msa.equals(new Object()));
  }

  public void testToString() {
    final MapSet<String, Integer> ms = new MapSet<>();
    ms.put("a", 1);
    ms.put("b", 2);
    ms.put("a", 1);
    assertEquals("{a=[1], b=[2]}", ms.toString());
  }
}

