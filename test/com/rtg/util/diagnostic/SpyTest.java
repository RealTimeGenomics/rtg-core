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
public class SpyTest extends TestCase {

  public void test() {
    Diagnostic.clearListeners();
    Diagnostic.setLogStream();
    final SpyCounter spy1 = new SpyCounter("foo1");
    final SpyCounter spy2 = new SpyCounter("bar2");

    spy1.increment();
    spy2.increment();
    spy2.increment();

    final MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());
    Spy.report();
    mps.close();
    final String t = mps.toString().trim();
    //System.err.println(t);
    assertTrue(t.contains("foo1 counts 1"));
    assertTrue(t.contains("bar2 counts 2"));
  }

}
