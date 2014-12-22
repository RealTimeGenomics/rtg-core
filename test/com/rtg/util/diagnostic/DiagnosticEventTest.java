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

import com.rtg.util.EnumHelper;
import com.rtg.util.PseudoEnum;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 * Tests the corresponding class.
 *
 */
public final class DiagnosticEventTest extends AbstractDiagnosticEventTest {

  /**
   */
  public DiagnosticEventTest(final String name) {
    super(name);
  }
  public static Test suite() {
    return new TestSuite(DiagnosticEventTest.class);
  }
  /**
   * Main to run from tests from command line.
   * @param args ignored.
   */
  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }

  @Override
  public DiagnosticEvent<?> getEvent() {
    return new DiagnosticEvent<ErrorType>(ErrorType.NOT_A_DIRECTORY, "1") { };
  }

  @Override
  public Class<?> getEnumClass() {
    return ErrorType.class;
  }

  public void testExplicitMessageForProgress() {
    assertTrue(getEvent().getMessage().contains("\"1\""));
  }

  public void test() {
    try {
      new DiagnosticEvent<>(ErrorType.NOT_A_DIRECTORY);
    } catch (IllegalArgumentException e) {
      assertEquals(ErrorType.NOT_A_DIRECTORY + ":" + 0 + ":" + 1, e.getMessage());
    }
  }

  static final class MyErrorType implements DiagnosticType, PseudoEnum {
    public static final MyErrorType NO_SUCH_ERROR = new MyErrorType(0, "NO_SUCH_ERROR");

    private final int mOrdinal;
    private final String mName;
    private MyErrorType(final int ordinal, final String name) {
      mOrdinal = ordinal;
      mName = name;
    }
    @Override
    public int ordinal() {
      return mOrdinal;
    }
    @Override
    public String name() {
      return mName;
    }

    public String toString() {
      return mName;
    }
    private static final EnumHelper<MyErrorType> HELPER = new EnumHelper<>(MyErrorType.class, new MyErrorType[] {NO_SUCH_ERROR});

    public static MyErrorType valueOf(final String str) {
      return HELPER.valueOf(str);
    }

    public static MyErrorType[] values() {
      return HELPER.values();
    }
    @Override
    public int getNumberOfParameters() {
      return 0;
    }
    @Override
    public String getMessagePrefix() {
      return "";
    }
  }

  public void testNonExistentErrorType() throws IOException {
    try {
      final ByteArrayOutputStream bos = new ByteArrayOutputStream();
      try {
        try (PrintStream ps = new PrintStream(bos)) {
          Diagnostic.setLogStream(ps);
          assertEquals("NO_SUCH_ERROR", new DiagnosticEvent<MyErrorType>(MyErrorType.NO_SUCH_ERROR) {
          }.getMessage());
        }
      } finally {
        bos.close();
      }
      final String t = bos.toString().trim();
      assertTrue(t.endsWith("Missing resource information for diagnostic: NO_SUCH_ERROR"));
      assertTrue(t.startsWith("20"));
      assertEquals('-', t.charAt(4));
      assertEquals('-', t.charAt(7));
      assertEquals(' ', t.charAt(10));
      assertEquals(':', t.charAt(13));
      assertEquals(':', t.charAt(16));
      assertEquals(' ', t.charAt(19));
    } finally {
      Diagnostic.setLogStream();
    }
  }
}

