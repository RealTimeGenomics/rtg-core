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
package com.rtg.simulation.reads;

import java.io.IOException;

import junit.framework.TestCase;

/**
 */
public class ErrorMachineTest extends TestCase {

  public void testDuplication() throws Exception {

    final int[] counts = new int[1];

    final Machine testMachine = new DummyMachineTest.MockMachine() {
      @Override
      public void processFragment(String id, int fragmentStart, byte[] data, int length) throws IOException {
        assertEquals((counts[0] > 0 ? "dupe-" /*+ (counts[0] + 1))*/ : "") + "blah1", id);
        counts[0]++;
      }

    };

    final ErrorMachine em = new ErrorMachine(0, testMachine, 1.0, 0.0);

    em.processFragment("blah1", 0, new byte[] {0, 1}, 1);

    assertEquals(2, counts[0]);
  }


  public void testChimera() throws Exception {

    final int[] counts = new int[1];

    final Machine testMachine = new DummyMachineTest.MockMachine() {
      @Override
      public void processFragment(String id, int fragmentStart, byte[] data, int length) throws IOException {
        assertEquals("chimera" + counts[0] + "/", id);
        assertEquals(Integer.MIN_VALUE, fragmentStart);
        assertEquals(3, length);
        assertEquals(3, data.length);
        assertEquals(0, data[0]);
        assertEquals(3, data[1]);
        assertEquals(2, data[2]);
        counts[0]++;
      }
    };

    final ErrorMachine em = new ErrorMachine(0, testMachine, 0.0, 1.0);

    em.processFragment("blah1", 0, new byte[] {0, 1}, 1);
    em.processFragment("blah2", 500, new byte[] {3, 2}, 2);

    em.processFragment("blah1", 0, new byte[] {0, 1}, 1);
    em.processFragment("blah2", 500, new byte[] {3, 2}, 2);

    assertEquals(2, counts[0]);
  }
}
