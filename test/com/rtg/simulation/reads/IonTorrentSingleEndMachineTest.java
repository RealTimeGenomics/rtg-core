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

import com.rtg.mode.DnaUtils;
import com.rtg.reader.SdfId;
import com.rtg.util.StringUtils;

import junit.framework.TestCase;

/**
 */
public class IonTorrentSingleEndMachineTest extends TestCase {

  public void testHomopolymerLength() {

    final byte[] data = DnaUtils.encodeString("GAAACCTCGTAAGTTTTTTTCTTTCTCTC");

    assertEquals(1, IonTorrentSingleEndMachine.homopolymerLength(0, data, data.length));
    assertEquals(3, IonTorrentSingleEndMachine.homopolymerLength(1, data, data.length));
    assertEquals(3, IonTorrentSingleEndMachine.homopolymerLength(2, data, data.length));
    assertEquals(3, IonTorrentSingleEndMachine.homopolymerLength(3, data, data.length));
    assertEquals(2, IonTorrentSingleEndMachine.homopolymerLength(4, data, data.length));
    assertEquals(2, IonTorrentSingleEndMachine.homopolymerLength(5, data, data.length));
    assertEquals(1, IonTorrentSingleEndMachine.homopolymerLength(6, data, data.length));
    assertEquals(1, IonTorrentSingleEndMachine.homopolymerLength(7, data, data.length));
    assertEquals(1, IonTorrentSingleEndMachine.homopolymerLength(8, data, data.length));
    assertEquals(1, IonTorrentSingleEndMachine.homopolymerLength(9, data, data.length));
    assertEquals(2, IonTorrentSingleEndMachine.homopolymerLength(10, data, data.length));
    assertEquals(2, IonTorrentSingleEndMachine.homopolymerLength(11, data, data.length));
    assertEquals(1, IonTorrentSingleEndMachine.homopolymerLength(12, data, data.length));
    assertEquals(7, IonTorrentSingleEndMachine.homopolymerLength(13, data, data.length));

  }


  public void testRandomHomopoly() throws Exception {
    final IonTorrentSingleEndMachine m = new IonTorrentSingleEndMachine(1312349873);
    m.reseedErrorRandom(1312349873);
    m.setMinSize(5);
    m.setMaxSize(20);
    m.updateWorkingSpace(10);

    final byte[] data = DnaUtils.encodeString("GAAACCTCGTAAGTTTTTTTCTTTCTCTC");
    int pos = m.readBases(0, data, data.length, 1, 10, 0, 1);

    assertEquals(11, pos);
    assertEquals("3.1D7.", m.getCigar(false));
    assertEquals("GAACCTCGTA", DnaUtils.bytesToSequenceIncCG(m.mReadBytes)); //reduced the A homopoly by 1

    String s = m.formatActionsHistogram();
    assertEquals("Total action count:\t11" + StringUtils.LS
        + "Match count:\t10\t90.91%" + StringUtils.LS
        + "Deletion count:\t1\t9.09%" + StringUtils.LS
        + "Total error count:\t1\t9.09%" + StringUtils.LS
        + "Of deletions, due to homopolymer:\t1\t9.09%" + StringUtils.LS, s);
    m.reseedErrorRandom(1348);
    m.resetCigar();
    pos = m.readBases(0, data, data.length, 1, 10, 0, 1);
    assertEquals("4.1I4.1I", m.getCigar(false));
    assertEquals("GAAACCCTCT", DnaUtils.bytesToSequenceIncCG(m.mReadBytes)); //elongated the A homopoly by 1, C by 1 (one of these is NOT due to homopolymer extension!)
    assertEquals(8, pos);
    s = m.formatActionsHistogram();
    assertEquals("Total action count:\t21" + StringUtils.LS
        + "Match count:\t18\t85.71%" + StringUtils.LS
        + "Deletion count:\t1\t4.76%" + StringUtils.LS
        + "Insertion count:\t2\t9.52%" + StringUtils.LS
        + "Total error count:\t3\t14.29%" + StringUtils.LS
        + "Of deletions, due to homopolymer:\t1\t4.76%" + StringUtils.LS
        + "Of insertions, due to homopolymer:\t1\t4.76%" + StringUtils.LS, s);
  }

  public void testEndN() throws Exception {
    final IonTorrentSingleEndMachine m = new IonTorrentSingleEndMachine(4);
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
    public void writeRead(String name, byte[] data, byte[] qual, int length) {
      mLastData = new byte[length];
      System.arraycopy(data, 0, mLastData, 0, length);
      mName = name;
    }

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
