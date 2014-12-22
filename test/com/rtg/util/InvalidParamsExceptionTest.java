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

import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.DiagnosticEvent;
import com.rtg.util.diagnostic.DiagnosticListener;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.io.LogRecord;
import com.rtg.util.io.LogStream;

import junit.framework.TestCase;

/**
 */
public class InvalidParamsExceptionTest extends TestCase {

  public void test() {
    final LogStream rec = new LogRecord();
    Diagnostic.setLogStream(rec);
    try {
      throw new InvalidParamsException(ErrorType.SLIM_ERROR);
    } catch (final InvalidParamsException e) {
      assertEquals(ErrorType.SLIM_ERROR, e.getErrorType());
    }
    //System.err.println(rec.toString());
  }

  public void testOtherConstructor() {
    final LogStream rec = new LogRecord();
    Diagnostic.setLogStream(rec);
    DiagnosticListener listener = new DiagnosticListener() {
      @Override
      public void close() {
      }
      @Override
      public void handleDiagnosticEvent(DiagnosticEvent<?> event) {
        assertEquals(ErrorType.INFO_ERROR, event.getType());
      }
    };
    Diagnostic.addListener(listener);
    try {
      try {
        throw new InvalidParamsException("This is a message");
      } catch (final InvalidParamsException e) {
        assertEquals("This is a message", e.getMessage());
      }
    } finally {
      Diagnostic.removeListener(listener);
    }
  }
}
