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

package com.rtg.util.diagnostic;

import com.rtg.util.io.MemoryPrintStream;

import junit.framework.TestCase;

/**
 */
public class SpyHistogramTest extends TestCase {

  public void test() {
    final MemoryPrintStream ps = new MemoryPrintStream();
    Diagnostic.setLogStream(ps.printStream());
    final SpyHistogram spy = new SpyHistogram("foo", 3);
    spy.increment(1);
    spy.increment(1);
    spy.increment(2);
    spy.increment(5);
    Spy.report();
    ps.close();
    //System.err.println(ps.toString());
    final String exp = "foo [3] 0 2 1 ...1";
    assertEquals(exp, spy.toString());
    assertTrue(ps.toString().contains(exp));
  }
}
