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
public class FastqReadWriterTest extends TestCase {

  public void test() throws Exception {
    final MemoryPrintStream out = new MemoryPrintStream();
    final FastqReadWriter f = new FastqReadWriter(out.printStream());
    f.writeRead("foo", new byte[] {0, 1, 2, 3, 4}, new byte[] {20, 20, 20, 20, 20}, 5);
    f.writeLeftRead("foo", new byte[] {0, 1, 2, 3, 4}, new byte[] {20, 20, 20, 20, 20}, 5);
    f.writeRightRead("foo", new byte[] {0, 1, 2, 3, 4}, new byte[] {20, 20, 20, 20, 20}, 5);
    assertEquals("@0 foo\nNACGT\n+\n55555\n@1 foo 1\nNACGT\n+\n55555\n@1 foo 2\nNACGT\n+\n55555\n", out.toString());
    assertEquals(2, f.readsWritten());
  }

}
