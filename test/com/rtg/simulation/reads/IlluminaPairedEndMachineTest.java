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

import com.rtg.util.io.MemoryPrintStream;

import junit.framework.TestCase;

/**
 * Test class
 */
public class IlluminaPairedEndMachineTest extends TestCase {

  public void testProcessFragment() throws Exception {
    final IlluminaPairedEndMachine m = new IlluminaPairedEndMachine(42);
    final MemoryPrintStream out = new MemoryPrintStream();
    final FastaReadWriter w = new FastaReadWriter(out.printStream());
    m.setReadWriter(w);
    assertTrue(m.isPaired());
    m.setLeftReadLength(5);
    m.setRightReadLength(5);
    final byte[] frag = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    m.processFragment("name/", 30, frag, frag.length);
    assertEquals(">0 name/31/F/5./Left\nAAAAA\n>0 name/52/R/5./Right\nTTTTT\n", out.toString());
  }

  public void testProcessFragmentReadThrough() throws Exception {
    final IlluminaPairedEndMachine m = new IlluminaPairedEndMachine(42);
    final MemoryPrintStream out = new MemoryPrintStream();
    final FastaReadWriter w = new FastaReadWriter(out.printStream());
    m.setReadWriter(w);
    m.setLeftReadLength(55);
    m.setRightReadLength(55);
    final byte[] frag = {1, 1, 1, 1, 1};
    m.processFragment("name/", 30, frag, frag.length);
    assertEquals(">0 name/31/F/55./Left\nAAAAACGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCTNNNNNTCTAGCCT\n>0 name/-19/R/55./Right\nTTTTTTGTGAGAAAGGGATGTGCTGCGAGAAGGCTAGANNNNNAGATCGGAAGAG\n", out.toString());
  }

  public void testBaseQuality() throws Exception {
    final IlluminaPairedEndMachine m = new IlluminaPairedEndMachine(42);
    m.setQualRange((byte) 15, (byte) 35);
    int qsummatch = 0;
    int qsummismatch = 0;
    final int n = 1000;
    for (int i = 0; i < n; i++) {
      final byte correctCallQuality = m.getCorrectCallQuality((byte) 1);
      final byte missCallQuality = m.getMissCallQuality();
      //System.err.println(String.format(" = %2d  X %2d  %b", correctCallQuality, missCallQuality, correctCallQuality > missCallQuality));
      qsummatch += correctCallQuality;
      qsummismatch += missCallQuality;
    }
//    System.err.println("Avg match quality = " + qsummatch / n);
//    System.err.println("Avg mismatch quality = " + qsummismatch / n);
    assertTrue(qsummatch > qsummismatch);
  }
}
