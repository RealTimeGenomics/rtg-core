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
package com.rtg.position;

import java.io.File;

import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.SlimException;
import com.rtg.util.io.LogRecord;
import com.rtg.util.io.LogStream;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 */
public class PositionUtilsTest extends TestCase {

  public static TestSuite suite() {
    final TestSuite suite = new TestSuite();
    suite.addTestSuite(com.rtg.position.PositionUtilsTest.class);
    return suite;
  }

  protected File mDir;

  @Override
  public void setUp() throws Exception {
    Diagnostic.setLogStream();
    mDir = FileHelper.createTempDirectory();
  }

  @Override
  public void tearDown() {
    assertTrue(FileHelper.deleteAll(mDir));
    mDir = null;
  }


  public void testMakeBuffer() {
    final byte[] ba = PositionUtils.makeBuffer(10);
    assertEquals(10, ba.length);
    final LogStream log = new LogRecord();
    Diagnostic.setLogStream(log);
    try {
      PositionUtils.makeBuffer(Integer.MAX_VALUE + 1L);
      fail();
    } catch (final SlimException e) {
      e.logException();
    } finally {
      Diagnostic.setLogStream();
    }
    final String str = log.toString();
    assertTrue(str.contains("There is a sequence which is too long to process. Its length is \"2147483648\" bytes. See the SDF output for the name of the sequence."));
    assertTrue(str.contains(SlimException.class.getName()));
    assertTrue(str.contains(" RTG has encountered a difficulty, please contact support@"));
  }

}




