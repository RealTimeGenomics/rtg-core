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
package com.rtg.util;

import com.rtg.util.ProgramState.SlimAbortException;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.LogRecord;

import junit.framework.TestCase;

/**
 */
public class ProgramStateTest extends TestCase {

  public void test() {
    final LogRecord log = new LogRecord();
    Diagnostic.setLogStream(log);
    //be careful this is all static stuff so cant be sure of the state coming in.
    final String msg = "Aborting operation in thread: ";
    assertFalse(log.toString().contains(msg));
    ProgramState.checkAbort();
    assertFalse(log.toString().contains(msg));
    ProgramState.setAbort();
    try {
      ProgramState.checkAbort();
      fail();
    } catch (final SlimAbortException e) {
      assertTrue(e.getMessage().contains(msg));
      assertTrue(log.toString().contains(msg));
    }
    Diagnostic.setLogStream();
    final LogRecord log2 = new LogRecord();
    Diagnostic.setLogStream(log2);
    ProgramState.clearAbort();
    ProgramState.checkAbort();
    assertFalse(log2.toString().contains(msg));
    Diagnostic.setLogStream();
  }
}
