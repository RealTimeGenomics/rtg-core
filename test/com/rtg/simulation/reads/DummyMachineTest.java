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

import com.rtg.reader.PrereadType;
import com.rtg.simulation.reads.AbstractMachine.SimErrorType;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.TestUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.variant.AbstractMachineErrorParams;
import com.rtg.variant.MachineErrorParamsBuilder;

import junit.framework.TestCase;

/**
 * Test class
 */
public class DummyMachineTest extends TestCase {

  public AbstractMachine getMachine() throws IOException, InvalidParamsException {
    return new AbstractMachine(new MachineErrorParamsBuilder().errors("illumina").create()) {
      @Override
      public void setReadWriter(ReadWriter rw) {
      }
      @Override
      public void processFragment(String id, int fragmentStart, byte[] data, int length) {
      }
      @Override
      public boolean isPaired() {
        return false;
      }
      @Override
      public PrereadType machineType() {
        return null;
      }
    };
  }

  public void test() throws Exception {
    final AbstractMachine m = getMachine();
    final MemoryPrintStream out = new MemoryPrintStream();
    final FastaReadWriter w = new FastaReadWriter(out.printStream());
    m.setReadWriter(w);
    m.mWorkspace = new int[10];
    m.mReadBytes = new byte[10];
    m.mQualityBytes = new byte[10];
    final byte[] frag = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    m.reseedErrorRandom(21);
    assertEquals(0, m.process(0, frag, frag.length, 1, 10));
    assertTrue(Arrays.equals(new byte[] {1, 1, 1, 1, 1, 1, 3, 1, 1, 1}, m.mReadBytes));

    assertEquals("6.1X3.", m.getCigar(false, 0, 10, 10));
    assertEquals("6.1X6.", m.getCigar(true, 0, 10, 10));

    assertEquals(0, m.mResidueCount);   //this is never incremented by process...
    assertEquals(0, m.residues());

    assertEquals(0, m.process(0, frag, frag.length, 1, 10));
    assertEquals("10.", m.getCigar(false, 0, 10, 10));
    assertEquals(0, m.process(0, frag, frag.length, 1, 10));
    assertEquals("1X9.", m.getCigar(false, 0, 10, 10));
    for (int c = 0; c < 1551; c++) {
      m.process(0, frag, frag.length, 1, 10);
    }
    //    System.out.println(m.getCigar(false, 0, 10, 10));
    assertEquals("9.1D1.", m.getCigar(false, 0, 10, 10));
    for (int c = 0; c < 1134; c++) {
      m.process(0, frag, frag.length, 1, 10);
    }
    assertEquals("3.1I6.", m.getCigar(false, 0, 10, 10));
  }

  public void testNameFormat() throws Exception {
    assertEquals("34/9/F/33.", AbstractMachine.formatReadName("34/", 'F', "33.", 5, 3));
    assertEquals("blah/F/33.", AbstractMachine.formatReadName("blah/", 'F', "33.", -1, 3));
  }

  public void testErrorType() {
    TestUtils.testEnum(SimErrorType.class, "[MNP, INSERT, DELETE, SUBSTITUTE_N, NOERROR]");

    final AbstractMachineErrorParams priors = new MachineErrorParamsBuilder().errorInsEventRate(0.2).errorDelEventRate(0.15).errorMnpEventRate(0.3).create();

    final AbstractMachine am = new AbstractMachine(priors) {
      @Override
      public void setReadWriter(ReadWriter rw) {
      }
      @Override
      public void processFragment(String id, int fragmentStart, byte[] data, int length) {
      }
      @Override
      public boolean isPaired() {
        return false;
      }
      @Override
      public PrereadType machineType() {
        return null;
      }
    };

    assertEquals(SimErrorType.INSERT, am.getErrorType(0.0));
    assertEquals(SimErrorType.INSERT, am.getErrorType(0.1999));
    assertEquals(SimErrorType.DELETE, am.getErrorType(0.20));
    assertEquals(SimErrorType.DELETE, am.getErrorType(0.3499));
    assertEquals(SimErrorType.MNP, am.getErrorType(0.35));
    assertEquals(SimErrorType.MNP, am.getErrorType(0.6499));
    assertEquals(SimErrorType.NOERROR, am.getErrorType(0.65));
    assertEquals(SimErrorType.NOERROR, am.getErrorType(0.99));

  }
}
