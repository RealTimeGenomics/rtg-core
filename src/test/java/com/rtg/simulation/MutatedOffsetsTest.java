/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
