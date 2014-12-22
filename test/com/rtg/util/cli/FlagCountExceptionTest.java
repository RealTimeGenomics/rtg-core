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
package com.rtg.util.cli;


import junit.framework.TestCase;

/**
 *
 *
 *
 */
public class FlagCountExceptionTest extends TestCase {

  /**
   */
  public FlagCountExceptionTest(final String name) {
    super(name);
  }

  /**
   * Test method for {@link com.rtg.util.cli.FlagCountException#FlagCountException(java.lang.String)}.
   */
  public final void testFlagCountException() {
    FlagCountException fex = new FlagCountException("This is flag count exception");
    try {
      throw fex;
    } catch (final FlagCountException ex) {
      assertEquals("This is flag count exception", ex.getMessage());
    }
  }

}
