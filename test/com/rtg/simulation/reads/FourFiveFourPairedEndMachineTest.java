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

import com.rtg.util.InvalidParamsException;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.variant.AbstractMachineErrorParams;
import com.rtg.variant.MachineErrorParamsBuilder;

/**
 * Test class
 */
public class FourFiveFourPairedEndMachineTest extends AbstractMachineTest {

  @Override
  protected AbstractMachineErrorParams getPriors() throws IOException, InvalidParamsException {
    return new MachineErrorParamsBuilder().errors("ls454_pe").create();
  }

  @Override
  protected Machine getMachine(final long seed) throws IOException, InvalidParamsException {
    final FourFiveFourPairedEndMachine m = new FourFiveFourPairedEndMachine(getPriors(), seed);
    m.setMinPairSize(300);
    m.setMaxPairSize(350);
    return m;
  }

  public void test() throws Exception {
    final FourFiveFourPairedEndMachine m = (FourFiveFourPairedEndMachine) getMachine(47);
    final MemoryPrintStream out = new MemoryPrintStream();
    final FastaReadWriter w = new FastaReadWriter(out.printStream());
    m.setReadWriter(w);
    m.setMinPairSize(3);
    m.setMaxPairSize(7);
    final byte[] frag = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    m.processFragment("name/", 30, frag, frag.length);
    m.processFragment("name/", 30, frag, frag.length);
    m.processFragment("name/", 30, frag, frag.length);

    //System.out.println(out.toString());
    mNano.check("454pe-results.fa", out.toString(), false);
  }


}
