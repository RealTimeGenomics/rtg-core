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

import com.rtg.mode.DnaUtils;
import com.rtg.reader.SdfId;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.variant.AbstractMachineErrorParams;
import com.rtg.variant.MachineErrorParamsBuilder;

/**
 * Test corresponding class.
 */
public class FourFiveFourSingleEndMachineTest extends AbstractMachineTest {

  @Override
  protected AbstractMachineErrorParams getPriors() throws IOException, InvalidParamsException {
    return new MachineErrorParamsBuilder().errors("ls454_pe").create();
  }

  @Override
  protected Machine getMachine(final long seed) throws IOException, InvalidParamsException {
    final FourFiveFourSingleEndMachine m = new FourFiveFourSingleEndMachine(getPriors(), seed);
    m.setMinSize(300);
    m.setMaxSize(350);
    return m;
  }

  public void test() throws Exception {
    final FourFiveFourSingleEndMachine m = (FourFiveFourSingleEndMachine) getMachine(47);
    final MemoryPrintStream out = new MemoryPrintStream();
    final FastaReadWriter w = new FastaReadWriter(out.printStream());
    m.setReadWriter(w);
    m.setMinSize(3);
    m.setMaxSize(7);
    final byte[] frag = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    m.processFragment("name/", 30, frag, 5);
    m.processFragment("name/", 30, frag, 5);
    m.processFragment("name/", 30, frag, 5);

    //System.out.println(out.toString());
    mNano.check("454se-results.fa", out.toString(), false);
  }

  public void testEndN() throws Exception {
    final FourFiveFourSingleEndMachine m = new FourFiveFourSingleEndMachine(4);
    m.reseedErrorRandom(4);
    m.setMinSize(50);
    m.setMaxSize(50);
    m.updateWorkingSpace(10);

    final MockReadWriter mrw = new MockReadWriter();
    m.mReadWriter = mrw;

    final byte[] t = DnaUtils.encodeString("TCAGCTATTGTTCACCTTTCTTCTATACTGTATGTATGTCTCAGCAAGCTTGTGTTTGTTTGGTGGTTGGCTCCTCTATCTGTGGATGCATCAACTCCATNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN");
    m.processFragment("id/", 0, t, 100);
    assertEquals("TGTGTTTGTTTGGTGGTTGGCTCCTCTATCTGTGGATGCATCAACTCCAT", DnaUtils.bytesToSequenceIncCG(mrw.mLastData));
    assertEquals("id/51/F/50.", mrw.mName);
  }

  private static class MockReadWriter implements ReadWriter {
    byte[] mLastData;
    String mName;

    @Override
    public void close() {
    }

    @Override
    public void identifyTemplateSet(SdfId... templateIds) {
    }

    @Override
    public void identifyOriginalReference(SdfId referenceId) {
    }

    @Override
    public void writeRead(String name, byte[] data, byte[] qual, int length) {
      mLastData = new byte[length];
      System.arraycopy(data, 0, mLastData, 0, length);
      mName = name;
    }

    @Override
    public void writeLeftRead(String name, byte[] data, byte[] qual, int length) {
      writeRead(name, data, qual, length);
    }

    @Override
    public void writeRightRead(String name, byte[] data, byte[] qual, int length) {
      writeRead(name, data, qual, length);
    }
    @Override
    public int readsWritten() {
      return 1;
    }
  }
}
