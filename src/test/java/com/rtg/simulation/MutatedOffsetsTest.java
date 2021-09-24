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

package com.rtg.simulation;

import com.rtg.util.diagnostic.SlimException;

import junit.framework.TestCase;

/**
 * Test Class
 */
public class MutatedOffsetsTest extends TestCase {

  public void testEmpty() {
    final MutatedOffsets offsets = new MutatedOffsets();
    offsets.freeze();
    assertEquals(100, offsets.getOffsetPosition(100));
  }

  public void testBasic() {
    final MutatedOffsets offsets = new MutatedOffsets();

    offsets.addInsertion(1506, 1);
    offsets.addDeletion(3503, 3);
    offsets.addDeletion(9777, 2);

    offsets.freeze();

    assertEquals(100, offsets.getOffsetPosition(100));
    assertEquals(1506, offsets.getOffsetPosition(1506));
    assertEquals(1507, offsets.getOffsetPosition(1507));
    assertEquals(1507, offsets.getOffsetPosition(1508));
    assertEquals(3503, offsets.getOffsetPosition(3504));
    assertEquals(3507, offsets.getOffsetPosition(3505));
    assertEquals(9777, offsets.getOffsetPosition(9775));
    assertEquals(9780, offsets.getOffsetPosition(9776));
    assertEquals(10004, offsets.getOffsetPosition(10000));
  }

  public void testMoreComplex() {
    //Set initial array size so growth is tested
    final MutatedOffsets offsets = new MutatedOffsets(5);

    offsets.addInsertion(1506, 1);
    offsets.addInsertion(3503, 3);
    offsets.addDeletion(3532, 4);

    offsets.freeze();

    assertEquals(100, offsets.getOffsetPosition(100));
    assertEquals(1506, offsets.getOffsetPosition(1506));
    assertEquals(1507, offsets.getOffsetPosition(1507));
    assertEquals(1507, offsets.getOffsetPosition(1508));
    assertEquals(3503, offsets.getOffsetPosition(3504));
    assertEquals(3503, offsets.getOffsetPosition(3505));
    assertEquals(3504, offsets.getOffsetPosition(3506));
    assertEquals(3504, offsets.getOffsetPosition(3507));
    assertEquals(3504, offsets.getOffsetPosition(3508));
    assertEquals(3505, offsets.getOffsetPosition(3509));
    assertEquals(3532, offsets.getOffsetPosition(3536));
    assertEquals(3537, offsets.getOffsetPosition(3537));
    assertEquals(4000, offsets.getOffsetPosition(4000));
  }

  public void testEvenInsertion() {
    final MutatedOffsets offsets = new MutatedOffsets();

    offsets.addInsertion(1506, 8);

    offsets.freeze();

    assertEquals(1505, offsets.getOffsetPosition(1505));
    assertEquals(1506, offsets.getOffsetPosition(1506));
    assertEquals(1506, offsets.getOffsetPosition(1507));
    assertEquals(1506, offsets.getOffsetPosition(1508));
    assertEquals(1506, offsets.getOffsetPosition(1509));
    assertEquals(1506, offsets.getOffsetPosition(1510));
    assertEquals(1507, offsets.getOffsetPosition(1511));
    assertEquals(1507, offsets.getOffsetPosition(1512));
    assertEquals(1507, offsets.getOffsetPosition(1513));
    assertEquals(1507, offsets.getOffsetPosition(1514));
    assertEquals(1507, offsets.getOffsetPosition(1515));
    assertEquals(1508, offsets.getOffsetPosition(1516));
  }

  public void testExceptions() {
    final MutatedOffsets offsets = new MutatedOffsets();
    try {
      offsets.getOffsetPosition(1);
      fail("Expected that needs to be frozen before lookup.");
    } catch (SlimException e) {
      assertEquals("Cannot use to get offsets until frozen.", e.getMessage());
    }
    try {
      offsets.addDeletion(2, 1);
      offsets.addDeletion(1, 1);
      fail("Expected that order is required");
    } catch (IllegalArgumentException e) {
      assertEquals("Cannot add a position earlier than existing positions.", e.getMessage());
    }
    try {
      offsets.addInsertion(20, 1);
      offsets.addInsertion(20, 1);
      fail("Expected that order is required");
    } catch (IllegalArgumentException e) {
      assertEquals("Cannot add a position earlier than existing positions.", e.getMessage());
    }
    offsets.freeze();
    try {
      offsets.addDeletion(1, 1);
      fail("Expected that cannot add if frozen.");
    } catch (SlimException e) {
      assertEquals("Cannot add offsets once frozen.", e.getMessage());
    }
    try {
      offsets.addInsertion(1, 1);
      fail("Expected that cannot add if frozen.");
    } catch (SlimException e) {
      assertEquals("Cannot add offsets once frozen.", e.getMessage());
    }
  }
}
