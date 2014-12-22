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
import java.util.Arrays;

import com.rtg.util.io.MemoryPrintStream;

/**
 * Test class
 */
public class IlluminaSingleEndMachineTest extends DummyIlluminaMachineTest {

  public void testPaired() throws Exception {
    final IlluminaSingleEndMachine m = new IlluminaSingleEndMachine(42);
    assertFalse(m.isPaired());
  }

  @Override
  public void test() throws Exception {
    final IlluminaSingleEndMachine m = (IlluminaSingleEndMachine) getMachine(42);
    final MemoryPrintStream out = new MemoryPrintStream();
    final FastaReadWriter w = new FastaReadWriter(out.printStream()) {
        @Override
        public void writeRead(String name, byte[] data, byte[] qual, int length) throws IOException {
          super.writeRead(name, data, qual, length);
          assertEquals("[20, 20, 20, 20]", Arrays.toString(qual));
        }
      };
    m.setReadWriter(w);
    final byte[] frag = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    m.setReadLength(4);

    m.processFragment("name/", 30, frag, frag.length);
    assertEquals(">0 name/31/F/4.\nAAAA\n", out.toString());

    assertEquals(4, m.mResidueCount);
    assertEquals(4, m.mWorkspace.length);
  }
}
