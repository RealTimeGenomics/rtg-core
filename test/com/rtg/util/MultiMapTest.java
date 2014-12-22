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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;



/**
 * Test class for MultiMap.
 */
public class MultiMapTest extends TestCase {
  /**
   * Constructor (needed for JUnit)
   * @param name A string which names the object.
   */
  public MultiMapTest(final String name) {
    super(name);
  }

  private static final String EXPECTED = ""
    + "{" + StringUtils.LS
    + " 10 -> [5, 6, 7, 10]" + StringUtils.LS
    + "}" + StringUtils.LS
    ;

  public void testDefaultCon() {
    MultiMap<Integer, String> mm = new MultiMap<>();
    checkMap(mm, EXPECTED);
  }

  public void testCon() {
    MultiMap<Integer, String> mm = new MultiMap<>(new HashMap<Integer, Collection<String>>()
        , new MultiMapFactory<String>() {
      @Override
      public Collection<String> createCollection() {
        return new ArrayList<>();
      }
    });
    checkMap(mm, EXPECTED);
  }

  public void checkMap(final MultiMap<Integer, String> mm, final String expected) {
    assertEquals(0, mm.size());
    assertEquals(0, mm.values().size());
    assertNull(mm.get(5));
    mm.put(3, "1");
    mm.put(3, "2");
    mm.put(3, "3");
    mm.put(7, "1");
    mm.put(9, "1");
    mm.put(7, "2");
    assertEquals(Arrays.asList("1", "2", "3"), mm.get(3));
    assertEquals(Arrays.asList("1", "2"), mm.get(7));
    assertEquals(3, mm.size());
    assertEquals(3, mm.values().size());
    assertEquals(2, mm.remove(7).size());
    assertEquals(2, mm.size());
    assertEquals(2, mm.values().size());
    mm.clear();
    assertEquals(null, mm.get(3));
    assertEquals(0, mm.size());
    assertEquals(0, mm.values().size());
    mm.put(10, "2");
    assertEquals(1, mm.size());
    assertEquals(1, mm.values().size());
    assertEquals(Arrays.asList("2"), mm.get(10));
    assertEquals(1, mm.set(10, Arrays.asList("5", "6", "7")).size());
    assertEquals(1, mm.size());
    assertEquals(Arrays.asList("5", "6", "7"), mm.get(10));
    mm.put(10, "10");
    assertEquals(1, mm.size());
    assertEquals(1, mm.values().size());
    assertEquals(Arrays.asList("5", "6", "7", "10"), mm.get(10));
    assertEquals(expected, mm.toString());
  }

  public static Test suite() {
    return new TestSuite(MultiMapTest.class);
  }


  /**
   * Main to run from tests from command line.
   * @param args ignored.
   */
  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }
}

