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


import junit.framework.Test;
import junit.framework.TestSuite;

/**
 * Tests the corresponding class.
 *
 */
public class InformationEventTest extends AbstractDiagnosticEventTest {

  /**
   */
  public InformationEventTest(final String name) {
    super(name);
  }
  public static Test suite() {
    return new TestSuite(InformationEventTest.class);
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
    return new InformationEvent(InformationType.PROCESSING_ITEM_N_OF_N, "", "blah", "1", "3");
  }

  @Override
  public Class<?> getEnumClass() {
    return InformationType.class;
  }

}

