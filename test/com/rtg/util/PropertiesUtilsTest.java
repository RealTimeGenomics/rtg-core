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

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Properties;

import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;

import junit.framework.TestCase;


/**
 * Tests for Property Utility class.
 */
public class PropertiesUtilsTest extends TestCase {

  public void testGetPropertiesResource() throws InvalidParamsException, IOException {
    Diagnostic.setLogStream();
    final Properties pr = PropertiesUtils.getPriorsResource("human", PropertiesUtils.PropertyType.PRIOR_PROPERTY);
    assertNotNull(pr);
    final ByteArrayOutputStream log = new ByteArrayOutputStream();
    final PrintStream ps = new PrintStream(log);
    Diagnostic.setLogStream(ps);
    try {
      PropertiesUtils.getPriorsResource("non-existant", PropertiesUtils.PropertyType.PRIOR_PROPERTY);
      fail();
    } catch (InvalidParamsException e) {
      assertEquals(ErrorType.INFO_ERROR, e.getErrorType());
      assertTrue(e.getMessage().contains("Invalid prior option \"non-existant"));
//      ps.flush();
//      TestUtils.containsAll(log.toString(), "non-existant", /*"as a properties file for priors"*/ "Invalid prior option");
    } finally {
      Diagnostic.setLogStream();
    }
  }

//  public void testGetIntegerProperty() {
//    int x = PropertiesUtils.getIntegerProperty("propertyutils.blah", 123);
//    assertEquals(123, x);
//
//    //%System.Environment.SetEnvironmentVariable("propertyutils.blah", "2345");
//    System.setProperty("propertyutils.blah", "2345");
//    x = PropertiesUtils.getIntegerProperty("propertyutils.blah", 123);
//    assertEquals(2345, x);
//
//    //%System.Environment.SetEnvironmentVariable("propertyutils.blah", "abc");
//    System.setProperty("propertyutils.blah", "abc");
//    try {
//      PropertiesUtils.getIntegerProperty("propertyutils.blah", 123);
//      fail("Expected an Exception");
//    } catch (RuntimeException e) {
//      //e.printStackTrace();
//    }
//  }
}
