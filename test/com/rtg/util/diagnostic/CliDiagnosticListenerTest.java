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
import java.io.PrintStream;

import junit.framework.TestCase;


/**
 * Tests the corresponding class.
 *
 */
public class CliDiagnosticListenerTest extends TestCase {

  public void test() {
    final ByteArrayOutputStream byteout = new ByteArrayOutputStream();
    final PrintStream output = new PrintStream(byteout);
    final CliDiagnosticListener l = new CliDiagnosticListener(output);
    l.handleDiagnosticEvent(null);
    l.handleDiagnosticEvent(new ErrorEvent(ErrorType.SLIM_ERROR));
    output.close();
    final String s = byteout.toString();
    assertTrue(s.contains("RTG has encountered a difficulty, please contact"));
  }

}

