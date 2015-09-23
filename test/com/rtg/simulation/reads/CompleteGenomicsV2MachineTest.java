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

import com.rtg.mode.DnaUtils;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.PortableRandom;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.test.RandomDna;
import com.rtg.variant.AbstractMachineErrorParams;

/**
 * Test Class
 */
public class CompleteGenomicsV2MachineTest extends CompleteGenomicsV1MachineTest {

  @Override
  protected Machine getMachine(final long seed) throws IOException, InvalidParamsException {
    return new CompleteGenomicsV2Machine(getPriors(), seed);
  }

  @Override
  public void test() throws IOException, InvalidParamsException {
    final CompleteGenomicsMachine m = (CompleteGenomicsMachine) getMachine(42);
    final MemoryPrintStream out = new MemoryPrintStream();
    final FastaReadWriter w = new FastaReadWriter(out.printStream());
    m.setReadWriter(w);
    final byte[] frag = new byte[500];
    Arrays.fill(frag, (byte) 1);
    m.processFragment("name/", 30, frag, frag.length);
    checkQualities(m.mQualityBytes);
    m.processFragment("name/", 30, frag, frag.length);
    checkQualities(m.mQualityBytes);
    m.processFragment("name/", 30, frag, frag.length);
    checkQualities(m.mQualityBytes);
    m.processFragment("name/", 30, frag, frag.length);
    checkQualities(m.mQualityBytes);
    mNano.check("cg-v2-results.fa", out.toString(), false);
  }

  public void test2() throws IOException, InvalidParamsException {
    final CompleteGenomicsMachine m = (CompleteGenomicsMachine) getMachine(10);
    final MemoryPrintStream out = new MemoryPrintStream();
    final FastaReadWriter w = new FastaReadWriter(out.printStream());
    m.setReadWriter(w);
    final String template = RandomDna.random(500, new PortableRandom(33));
    final byte[] frag = DnaUtils.encodeString(template);
    m.processFragment("name/", 30, frag, frag.length);
    checkQualities(m.mQualityBytes);
    m.processFragment("name/", 30, frag, frag.length);
    checkQualities(m.mQualityBytes);
    m.processFragment("name/", 30, frag, frag.length);
    checkQualities(m.mQualityBytes);
    m.processFragment("name/", 30, frag, frag.length);
    checkQualities(m.mQualityBytes);
    m.processFragment("name/", 30, frag, frag.length);
    checkQualities(m.mQualityBytes);
    m.processFragment("name/", 30, frag, frag.length);
    checkQualities(m.mQualityBytes);
    mNano.check("cg-v2-results2.fa", out.toString(), false);
  }

  class StatsReadWriter extends CompleteGenomicsV1MachineTest.StatsReadWriter {
    private final int[] mBackstep = new int[10];

    @Override
    void augment(final String name) {
      final int b = name.indexOf('B');
      if (b == -1) {
        mBackstep[0]++;
      } else {
        mBackstep[name.charAt(b - 1) - '0']++;
      }
      final int n1 = name.indexOf('N');
      assertEquals(-1, n1);
      mTotal++;
    }

    @Override
    void performStatisticalTests() throws IOException, InvalidParamsException {
      assertEquals(0, mBackstep[8]);
      final AbstractMachineErrorParams errors = getPriors();
      checkDiscreteDistribution("Overlap", errors.overlapDistribution2(), new int[] {mBackstep[7], mBackstep[6], mBackstep[5], mBackstep[4], mBackstep[3], mBackstep[2], mBackstep[1], mBackstep[0]}, 1);
    }
  }

  @Override
  public void testOverlapDistributions() throws Exception {
    try (StatsReadWriter w = new StatsReadWriter()) {
      final Machine m = getMachine(System.currentTimeMillis());
      m.setReadWriter(w);
      final byte[] frag = new byte[FRAGMENT_LENGTH];
      Arrays.fill(frag, (byte) 1);
      for (int k = 0; k < 10000; k++) {
        m.processFragment("b/", 0, frag, frag.length);
      }
      w.performStatisticalTests();
    }
  }

}
