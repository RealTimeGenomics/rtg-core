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
package com.rtg.ngs.blocking;

import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.MemoryPrintStream;

/**
 * Test class
 */
public class ReadBlockerSyncTest extends ReadBlockerTest {

  @Override
  ReadBlocker getReadBlocker(int reads, int threshold) {
    return new ReadBlockerSync(reads, threshold, "test-sync");
  }

  @Override
  protected String expectedPairingsString() {
    return "test-sync";
  }

  @Override
  public void test() {
    final ReadBlocker b = getReadBlocker(1, 255);
    b.increment(0);
    assertEquals(1, b.getCount(0));
    b.reset(0);
    assertEquals(0, b.getCount(0));
    for (int i = 0; i < 65535; ++i) {
      b.increment(0);
    }
    final MemoryPrintStream ps = new MemoryPrintStream();
    Diagnostic.setLogStream(ps.printStream());
    b.close();
    assertTrue(ps.toString(), ps.toString().contains("1 reads had count 255"));
    Diagnostic.setLogStream();
    ps.close();
  }
}
