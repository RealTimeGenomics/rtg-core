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

import com.rtg.util.TestUtils;
import com.rtg.util.io.MemoryPrintStream;

import junit.framework.TestCase;

/**
 */
public class SpyTimerTest extends TestCase {

  public void test() {
    final MemoryPrintStream ps = new MemoryPrintStream();
    Diagnostic.setLogStream(ps.printStream());
    final SpyTimer spy = new SpyTimer("foo");
    spy.start();
    spy.stop();
    spy.start();
    spy.stop();
    Spy.report();
    ps.close();
    //System.err.println(ps.toString());
    final String str = spy.toString();
    TestUtils.containsAll(str, "Timer foo", "count 2");
    assertTrue(ps.toString().contains(str));
  }
}
