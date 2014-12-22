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

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;

import junit.framework.TestCase;

/**
 */
public class OneShotTimerTest extends TestCase {

  /**
   * Test method for {@link com.rtg.util.diagnostic.OneShotTimer#stopLog()}.
   * @throws IOException if an I/O error occurs.
   */
  public final void testStopLog() throws IOException {
    final ByteArrayOutputStream ba = new ByteArrayOutputStream();
    final PrintStream pr = new PrintStream(ba);
    Diagnostic.setLogStream(pr);
    final long start = System.nanoTime();
    final OneShotTimer ti = new OneShotTimer("One_timer");
    ti.stopLog();
    final long finish = System.nanoTime();

    final String s = ba.toString();
    //System.err.println(s);
    assertTrue(s.contains(" Timer One_timer "));
    double time = Double.parseDouble(s.substring(s.lastIndexOf(' ')));
    double timeActual = (finish - start) / 1000000000.0;
    assertEquals(time, timeActual, 0.1);
    pr.close();
    ba.close();
    Diagnostic.setLogStream();
  }

  public final void testFormat() {
    final OneShotTimer ti = new OneShotTimer("One_timer");
    assertEquals("Timer One_timer      0.12", ti.toString(123456789L));
    assertEquals("One_timer", ti.toString());
  }
}

