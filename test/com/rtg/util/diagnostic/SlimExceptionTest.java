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

import junit.framework.Assert;
import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;


/**
 * Tests the corresponding class.
 *
 */
public class SlimExceptionTest extends TestCase {

  //    public SlimExceptionTest(final String name) {
  //        super(name);
  //    }

  public static Test suite() {
    return new TestSuite(SlimExceptionTest.class);
  }

  /**
   * Main to run from tests from command line.
   * @param args ignored.
   */
  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }

  public void test() throws IOException {
    final DiagnosticListener l = new DiagnosticListener() {
      @Override
      public void handleDiagnosticEvent(final DiagnosticEvent<?> event) {
        Assert.assertTrue(event instanceof ErrorEvent);
        bump();
      }

      @Override
      public void close() {
      }
    };
    Diagnostic.addListener(l);
    checkExceptionLogging(new SlimException());
    assertEquals(1, mCount);
    checkExceptionLogging(new SlimException(new RuntimeException()));
    assertEquals(2, mCount);
    checkExceptionLogging(new SlimException(ErrorType.NOT_A_DIRECTORY, "dir"));
    assertEquals(3, mCount);
    checkExceptionLogging(new SlimException(new IOException(), ErrorType.NOT_A_DIRECTORY, "dir"));
    assertEquals(4, mCount);
    checkExceptionLogging(new SlimException(new OutOfMemoryError(), ErrorType.SLIM_ERROR));
    assertEquals(5, mCount);

    SlimException se = new SlimException(false, null, null);
    assertEquals("", se.getMessage());
    checkExceptionLogging(se);
    assertEquals(5, mCount);

    se = new SlimException(new RuntimeException(), ErrorType.INFO_ERROR, "blkjsdfk");
    final String expected = "java.lang.RuntimeException";
    assertTrue(se.getMessage().startsWith(expected));
    assertEquals(ErrorType.INFO_ERROR, se.getErrorType());
    checkExceptionLogging(se);
    assertEquals(6, mCount);

    se = new SlimException("bjseh!");
    assertEquals(ErrorType.INFO_ERROR, se.getErrorType());
    checkExceptionLogging(se);
    assertEquals(7, mCount);

    Diagnostic.removeListener(l);
  }

  private void checkExceptionLogging(final SlimException e) throws IOException {
    ByteArrayOutputStream bos = new ByteArrayOutputStream();
    try {
      final PrintStream ps = new PrintStream(bos);
      Diagnostic.setLogStream(ps);
      try {
        e.logException();
        e.printErrorNoLog();
      } finally {
        ps.close();
      }
    } finally {
      Diagnostic.setLogStream();
      bos.close();
    }
    checkLog(bos.toString().trim());
  }

  private int mCount = 0;

  private void bump() {
    mCount++;
  }

  private void checkLog(final String t) {
    assertTrue(t.length() > 0);
  }
}

