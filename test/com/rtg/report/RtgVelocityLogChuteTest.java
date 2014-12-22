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
package com.rtg.report;

import org.apache.velocity.runtime.log.LogChute;

import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.MemoryPrintStream;

import junit.framework.TestCase;

/**
 */
public class RtgVelocityLogChuteTest extends TestCase {

  public void test() {
    final MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());
    try {
      final RtgVelocityLogChute chute = new RtgVelocityLogChute();
      chute.log(LogChute.TRACE_ID, "not gonna log");
      chute.log(LogChute.DEBUG_ID, "totally logged");
      assertFalse(chute.isLevelEnabled(LogChute.TRACE_ID));
      assertTrue(chute.isLevelEnabled(LogChute.DEBUG_ID));
      assertTrue(chute.isLevelEnabled(LogChute.INFO_ID));
    } finally {
      Diagnostic.setLogStream();
    }
    assertFalse(mps.toString().contains("not gonna log"));
    assertTrue(mps.toString().contains("totally logged"));
  }
}
