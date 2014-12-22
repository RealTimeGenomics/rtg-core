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


import junit.framework.TestCase;

/**
 * Tests the corresponding class.
 *
 */
public abstract class AbstractDiagnosticEventTest extends TestCase {

  /**
   */
  public AbstractDiagnosticEventTest(final String name) {
    super(name);
  }

  public abstract DiagnosticEvent<?> getEvent();


  public abstract Class<?> getEnumClass();

  public void testEventParameters() {
    final DiagnosticEvent<?> event = getEvent();
    assertNotNull(event.getType());
    final String[] params = event.getParams();
    assertNotNull(params);
    if (params.length > 0) {
      params[0] = "hi_there_bilbo";
      assertFalse(params[0].equals(event.getParams()[0]));
    }
  }

  public void testMessageGeneration() throws Exception {
    final String methodName = "values";
    for (final DiagnosticType type : (DiagnosticType[]) getEnumClass().getMethod(methodName, new Class<?>[] {}).invoke(null)) {
      final String[] params = new String[type.getNumberOfParameters()];
      for (int k = 0; k < params.length; k++) {
        params[k] = "XX" + k + "XX";
      }
      final String message = new DiagnosticEvent<>(type, params).getMessage();
      assertNotNull(message);
      assertTrue(message.length() > 0);
      for (String param : params) {
        assertTrue(message, message.contains(param));
      }
    }
  }

}

