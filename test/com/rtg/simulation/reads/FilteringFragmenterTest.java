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
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SdfId;
import com.rtg.reader.SequencesReader;
import com.rtg.simulation.genome.SequenceDistribution;
import com.rtg.util.StringUtils;
import com.rtg.util.intervals.ReferenceRegions;

import junit.framework.TestCase;

/**
 */
public class FilteringFragmenterTest extends TestCase {
  private static class MockMachine implements Machine {
    Integer mLastStart;
    @Override
    public void setQualRange(byte minq, byte maxq) {
    }

    @Override
    public void setReadWriter(ReadWriter rw) {
    }

    @Override
    public void identifyTemplateSet(SdfId... templateIds) {
    }

    @Override
    public void identifyOriginalReference(SdfId referenceId) {
    }

    @Override
    public void processFragment(String id, int fragmentStart, byte[] data, int length) {
      mLastStart = fragmentStart;
    }

    @Override
    public long residues() {
      return 0;
    }

    @Override
    public boolean isPaired() {
      return false;
    }

    @Override
    public PrereadType machineType() {
      return null;
    }

    @Override
    public String formatActionsHistogram() {
      return null;
    }
  }
  public void test() throws IOException {
    final ReferenceRegions bed = new ReferenceRegions();
    bed.add("foo", 20, 29);
    final SequenceDistribution[] dist = {new SequenceDistribution(new double[] {1.0})};
    final SequencesReader[] readers = {ReaderTestUtils.getReaderDnaMemory(">foo" + StringUtils.LS + "ACGTACCCACAGAGATAGACACACGTAGATGACACAGCCATGTCCCGCCATAT")};
    final MockMachine m = new MockMachine();
    final FilteringFragmenter fragmenter = new FilteringFragmenter(bed, 23, dist, readers);
    fragmenter.setMinFragmentSize(2);
    fragmenter.setMaxFragmentSize(2);
    fragmenter.setMachine(m);
    fragmenter.emitFragment(2, 0, 0, "foo", 18);
    assertNull(m.mLastStart);
    fragmenter.emitFragment(2, 0, 0, "foo", 19);
    assertEquals(Integer.valueOf(19), m.mLastStart);
    fragmenter.emitFragment(2, 0, 0, "foo", 29);
    assertEquals(Integer.valueOf(19), m.mLastStart);
    fragmenter.emitFragment(2, 0, 0, "foo", 28);
    assertEquals(Integer.valueOf(28), m.mLastStart);
  }
}
