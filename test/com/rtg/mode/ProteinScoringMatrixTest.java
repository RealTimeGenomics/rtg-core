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
package com.rtg.mode;

import java.io.IOException;
import java.util.MissingResourceException;

import com.rtg.util.InvalidParamsException;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.io.LogRecord;
import com.rtg.util.io.LogStream;

import junit.framework.TestCase;

/**
 */
public class ProteinScoringMatrixTest extends TestCase {

  private ProteinScoringMatrix mMatrix;

  /**
   * Set up the matrix
   */
  @Override
  public void setUp() throws IOException, InvalidParamsException {
    mMatrix = new ProteinScoringMatrix();
  }

  /**
   * Tear down the matrix
   */
  @Override
  public void tearDown() {
    mMatrix = null;
  }

  /**
   * Test method for {@link Blosum62#score(int, int, String)}.
   */
  public final void testScore() {
    check(4, Protein.A, Protein.A);
    check(-1, Protein.A, Protein.R);
    check(1, Protein.STOP, Protein.STOP);
    check(-1, Protein.K, Protein.D);
    check(2, Protein.M, Protein.L);
    check(0, Protein.X, Protein.A);
    check(-1, Protein.X, Protein.R);
    check(-2, Protein.X, Protein.C);
    check(-4, Protein.X, Protein.STOP);
  }

  private void check(final int sc, final Protein a, final Protein b) {
    assertEquals(sc, mMatrix.score(a.ordinal(), b.ordinal()));
    assertEquals(sc, mMatrix.score(b.ordinal(), a.ordinal()));
  }

  public final void testToStringStringBuilder() {
    final String str = mMatrix.toString();
    //System.err.println(str);
    assertTrue(str.startsWith("[-1, -4, 0, -1, -1, -1, -2,"));
    assertTrue(str.endsWith("-1, -2, -2, 0, -3, -1, 4]" + StringUtils.LS));
  }

  public final void testGetters() {
    assertEquals("K", 0.0410, mMatrix.getK());
    assertEquals("EXPECTED", -0.5209, mMatrix.getExpected());
    assertEquals("LAMBDA", 0.267, mMatrix.getLambda());
    assertEquals("HIT", 4.5, mMatrix.getHit());
    assertEquals("MISS", -1.0, mMatrix.getMiss());
    assertEquals("GAP", -10.0, mMatrix.getGap());
    assertEquals("EXTEND", -1.0, mMatrix.getExtend());
    assertEquals(11, mMatrix.getMaxScore());
  }

  public final void testConstructor() throws InvalidParamsException, IOException {
    try {
      new ProteinScoringMatrix("NOTEXISTS");
      fail();
    } catch (final MissingResourceException e) {
      assertTrue(e.getMessage(), e.getMessage().startsWith("Could not find:com/rtg/mode/NOTEXISTS"));
    }
  }

  public final void testConstructor2() throws InvalidParamsException, IOException {
    try {
      new ProteinScoringMatrix("BLOSUM62CORRUPT");
      fail();
    } catch (final MissingResourceException e) {
      assertTrue(e.getMessage(), e.getMessage().startsWith("Malformed resource: com/rtg/mode/BLOSUM62CORRUPT message:"));
    }
  }

  public final void testConstructor3() throws InvalidParamsException, IOException {
    try {
      new ProteinScoringMatrix("BLOSUM45TESTPR");
      fail();
    } catch (final MissingResourceException e) {
      assertTrue(e.getMessage(), e.getMessage().startsWith("Could not find:com/rtg/mode/BLOSUM45TESTPR.properties"));
    }
  }

  public final void testConstructor4() throws InvalidParamsException, IOException {
    final LogStream logStream = new LogRecord();
    Diagnostic.setLogStream(logStream);
    try {
      new ProteinScoringMatrix("BLOSUM45TEST");
      fail();
    } catch (final InvalidParamsException e) {
      assertEquals(ErrorType.PROPS_KEY_NOT_FOUND, e.getErrorType());
//      final String str = logStream.toString();
      //System.err.println("LOG:" + str);
//      TestUtils.containsAll(str, new String[] {"Unable to locate key"});
    } finally {
      Diagnostic.setLogStream();
    }
  }

  /**
   * @throws InvalidParamsException
   * @throws IOException
   */
  public final void testConstructor5() throws InvalidParamsException, IOException {
    final LogStream logStream = new LogRecord();
    Diagnostic.setLogStream(logStream);
    try {
      new ProteinScoringMatrix("BLOSUM45TESTUNI");
      fail();
    } catch (final InvalidParamsException e) {
      assertEquals(ErrorType.PROPS_INVALID, e.getErrorType());

//      final String str = logStream.toString();
//      final String exp = "Matrix property file \"com/rtg/mode/BLOSUM45TESTUNI.properties\" is invalid (contains illegal Unicode escape characters).";
//      assertTrue("Exception:" + e.getClass().getName() + " Messasge: " + e.getMessage() + " Actual: " + str + "\n" + "Expected to contain: " + exp, str.contains(exp));
    } finally {
      Diagnostic.setLogStream();
    }
  }
}
