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

package com.rtg.assembler.graph.implementation;

import com.rtg.assembler.graph.PathsIterator;

import junit.framework.TestCase;

/**
 */
public class IteratorLocalTest extends TestCase {

  public void test() {
    final GraphImplementation graph = GraphImplementationTest.graph();
    final IteratorLocal it = new IteratorLocal(graph);
    it.integrity();
    assertTrue(it.graph() == graph);
    checkState(it);
    assertEquals(0, it.nextPathId());
    checkState(it);

    it.set(+3);
    it.integrity();
    checkState(it);
    assertEquals(+4, it.nextPathId());
    assertEquals(0, it.contigIndex());
    assertEquals(-3, it.nextPathId());
    assertEquals(2, it.contigIndex());
    assertEquals(-2, it.nextPathId());
    assertEquals(1, it.contigIndex());
    assertEquals(0, it.nextPathId());
    checkState(it);
  }

  private void checkState(final PathsIterator it) {
    try {
      it.contigIndex();
      fail();
    } catch (final IllegalStateException e) {
      // TODO: handle exception
    }
  }
}
