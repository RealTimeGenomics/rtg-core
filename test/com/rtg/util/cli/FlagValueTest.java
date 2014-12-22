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
package com.rtg.util.cli;


import junit.framework.TestCase;

/**
 *
 *
 *
 */
public class FlagValueTest extends TestCase {

  /**
   */
  public FlagValueTest(final String name) {
    super(name);
  }

  public void testFlagIntVlaue() {
    Flag f = new Flag('I', "integer", "integer Flag", 1, 1, Integer.class, "Integer Value", Integer.MAX_VALUE, "");
    FlagValue v = new FlagValue(f, 55);
    assertEquals(55, v.getValue());
    assertEquals("integer=55", v.toString());
    Flag n = v.getFlag();
    assertEquals(f.isSet(), n.isSet());
    assertEquals(f.getCount(), n.getCount());
    assertEquals(f.getDescription(), n.getDescription());
    assertEquals(f.getMaxCount(), n.getMaxCount());
    assertEquals(f.getMinCount(), n.getMinCount());
    assertEquals(f.getChar(), n.getChar());
    assertEquals(f.getFlagUsage(), n.getFlagUsage());
    assertEquals(f.getCompactFlagUsage(), n.getCompactFlagUsage());
    assertTrue(f.equals(n));
  }

  public void testFlagDoubleVlaue() {
    Flag f = new Flag('D', "double", "double Flag", 1, 1, Double.class, "Double Value", Double.MAX_VALUE, "");
    FlagValue v = new FlagValue(f, 55.000256);
    assertEquals(55.000256, v.getValue());
    assertEquals("double=55.000256", v.toString());
    Flag n = v.getFlag();
    assertEquals(f.isSet(), n.isSet());
    assertEquals(f.getCount(), n.getCount());
    assertEquals(f.getDescription(), n.getDescription());
    assertEquals(f.getMaxCount(), n.getMaxCount());
    assertEquals(f.getMinCount(), n.getMinCount());
    assertEquals(f.getChar(), n.getChar());
    assertEquals(f.getFlagUsage(), n.getFlagUsage());
    assertEquals(f.getCompactFlagUsage(), n.getCompactFlagUsage());
    assertTrue(f.equals(n));
  }
}
