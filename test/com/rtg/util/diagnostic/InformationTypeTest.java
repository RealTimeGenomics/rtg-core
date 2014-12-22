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

import java.io.ObjectStreamException;

import com.rtg.util.TestUtils;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * Tests for InformationType.
 *
 */
public class InformationTypeTest extends TestCase {

  /**
   */
  public InformationTypeTest(final String name) {
    super(name);
  }
  public static Test suite() {
    return new TestSuite(InformationTypeTest.class);
  }
  /**
   * Main to run from tests from command line.
   * @param args ignored.
   */
  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }

  public void testEnum() {
    TestUtils.testPseudoEnum(InformationType.class, "[INFO_USER, PROCESSING_ITEM_N_OF_N]");
    assertEquals(1, InformationType.INFO_USER.getNumberOfParameters());
    assertEquals(4, InformationType.PROCESSING_ITEM_N_OF_N.getNumberOfParameters());
  }

  public void testReadResolve() throws ObjectStreamException {
    assertEquals(InformationType.INFO_USER, InformationType.INFO_USER.readResolve());
  }

  public void testPrefix() {
    assertEquals("", InformationType.INFO_USER.getMessagePrefix());
  }
}

