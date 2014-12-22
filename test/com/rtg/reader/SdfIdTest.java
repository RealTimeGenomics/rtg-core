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

package com.rtg.reader;

import java.util.UUID;

import junit.framework.TestCase;

/**
 * Test class
 */
public class SdfIdTest extends TestCase {

  public void testSomeMethod() {
    SdfId id = new SdfId();
    assertTrue(id.available());
    assertFalse(id.check(new SdfId()));
    assertTrue(id.check(new SdfId(0L)));
    assertFalse(id.check(new SdfId(1L)));
    assertFalse(new SdfId(0L).equals(id));
    assertEquals(new SdfId(new UUID(id.getHighBits(), id.getLowBits())), id);
    assertTrue(id.toString().matches("^[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{12}$"));
    assertTrue(0L != id.getLowBits());
    assertTrue(0L != id.getHighBits());
    id = new SdfId(555L);
    assertTrue(id.available());
    assertFalse(id.check(new SdfId()));
    assertTrue(id.check(new SdfId(0L)));
    assertFalse(id.check(new SdfId(1L)));
    assertFalse(new SdfId(0L).equals(id));
    assertFalse(new SdfId().equals(id));
    assertEquals(new SdfId(555L), id);
    assertTrue(id.toString().matches("^[1-9a-fA-F]([1-9a-fA-F]?){15}$"));
    assertEquals(555L, id.getLowBits());
    assertEquals(0L, id.getHighBits());
    id = new SdfId(0L);
    assertFalse(id.available());
    assertTrue(id.check(new SdfId()));
    assertTrue(id.check(new SdfId(0L)));
    assertTrue(id.check(new SdfId(1L)));
    assertEquals(new SdfId(0L), id);
    assertTrue(id.toString().matches("^0$"));
    assertEquals(0L, id.getLowBits());
    assertEquals(0L, id.getHighBits());
    id = new SdfId("ffffffff-ffff-ffff-ffff-ffffffffffff");
    assertFalse(id.check(new SdfId()));
    assertTrue(id.check(new SdfId(0L)));
    assertFalse(id.check(new SdfId(1L)));
    assertEquals("ffffffff-ffff-ffff-ffff-ffffffffffff", id.toString());
    assertEquals(-1L, id.getLowBits());
    assertEquals(-1L, id.getHighBits());
    id = new SdfId("ffffffff");
    assertFalse(id.check(new SdfId()));
    assertTrue(id.check(new SdfId(0L)));
    assertFalse(id.check(new SdfId(1L)));
    assertEquals("ffffffff", id.toString());
    assertEquals(0xffffffffL, id.getLowBits());
    assertEquals(0L, id.getHighBits());
    assertFalse(id.equals("1"));
  }

}
