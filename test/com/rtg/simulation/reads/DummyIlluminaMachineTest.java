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

import com.rtg.reader.PrereadType;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.variant.AbstractMachineErrorParams;
import com.rtg.variant.MachineErrorParamsBuilder;

/**
 * Test class
 */
public class DummyIlluminaMachineTest extends AbstractMachineTest {

  @Override
  protected AbstractMachineErrorParams getPriors() throws IOException, InvalidParamsException {
    return new MachineErrorParamsBuilder().errors("illumina").create();
  }

  @Override
  protected Machine getMachine(final long seed) throws IOException, InvalidParamsException {
    final IlluminaSingleEndMachine m = new IlluminaSingleEndMachine(getPriors(), seed);
    m.setReadLength(36);
    return m;
  }

  public void test() throws Exception {
    final AbstractIlluminaMachine m = (IlluminaSingleEndMachine) getMachine(42);
    final MemoryPrintStream out = new MemoryPrintStream();
    final FastaReadWriter w = new FastaReadWriter(out.printStream());
    m.setReadWriter(w);
    final byte[] frag = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    m.mWorkspace = new int[10];
    m.mReadBytes = new byte[10];
    m.reseedErrorRandom(m.mFrameRandom.nextLong());
    final String resF = m.generateRead("name/", 30, frag, frag.length, true, 5);
    assertEquals("name/31/F/5.", resF);
    final String resR = m.generateRead("name/", 30, frag, frag.length, false, 6);
    assertEquals("name/51/R/6.", resR);
    assertEquals(PrereadType.UNKNOWN, m.machineType());
  }
}
