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
package com.rtg.launcher;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.io.StringWriter;

import com.rtg.util.TestUtils;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.cli.TestCFlags;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.test.NanoRegression;

import junit.framework.TestCase;


/**
 * Abstract tests for AbstractCli subclasses
 *
 */
public abstract class AbstractCliTest extends TestCase {

  protected AbstractCli mCli;

  protected NanoRegression mNano;

  @Override
  public void setUp() throws IOException {
    GlobalFlags.resetAccessedStatus();
    CommandLine.clearCommandArgs();
    Diagnostic.setLogStream();
    mCli = getCli();
    mNano = new NanoRegression(this.getClass());
  }

  @Override
  public void tearDown() throws IOException {
    Diagnostic.setLogStream();
    mCli = null;
    try {
      mNano.finish();
    } finally {
      mNano = null;
    }
  }

  protected abstract AbstractCli getCli();

  protected CFlags getCFlags() {
    return mCli.mFlags;
  }

  public void testApplicationName() {
    assertEquals("rtg", mCli.applicationName());
  }

  /**
   * Checks the help output of the CLI class for consistency and that
   * it contains the specified strings.
   *
   * @param expected a <code>String</code> value
   */
  protected void checkHelp(String... expected) {
    mCli.createRegisterFlags(TestUtils.getNullPrintStream(), null);
    TestCFlags.check(mCli.getCFlags(), expected);
  }

  protected void checkExtendedHelp(String... expected) {
    mCli.createRegisterFlags(TestUtils.getNullPrintStream(), null);
    TestCFlags.checkExtendedUsage(mCli.getCFlags(), expected);
  }


  /**
   * Runs the supplied arguments through the CFlags, under the
   * assumption that they should fail validation.
   *
   * @param args command line arguments.
   * @return the concatenation of stderr and the log stream.
   */
  protected String checkHandleFlagsErr(String... args) {
    final ByteArrayOutputStream os = new ByteArrayOutputStream();
    final PrintStream ps = new PrintStream(os);
    Diagnostic.setLogStream(ps);
    final StringWriter err = new StringWriter();
    try {
      assertFalse(mCli.handleFlags(args, TestUtils.getNullPrintStream(), err));
      assertNotNull(mCli.mFlags);
      ps.flush();
    } finally {
      Diagnostic.setLogStream();
    }
    return os.toString() + err.toString();
  }

  /**
   * Runs the supplied arguments through the CFlags, under the
   * assumption that they should be validated fine.
   *
   * @param args command line arguments.
   * @return the contents of stdout.
   */
  protected String checkHandleFlagsOut(String... args) {
    final StringWriter writer = new StringWriter();
    final MemoryPrintStream err = new MemoryPrintStream();
    boolean val = mCli.handleFlags(args, writer, err.printStream());
    assertTrue(err.toString(), val);
    return writer.toString();
  }

  /**
   * Runs the supplied arguments through the CFlags.
   * @param args command line arguments.
   */
  protected void checkHandleFlags(String... args) {
    mCli.handleFlags(args, TestUtils.getNullPrintStream(), TestUtils.getNullPrintStream());
  }

  /**
   * Runs the main method of the cli class under the assumption that
   * it should run without errors.
   *
   * @param args command line arguments.
   * @return the contents of stdout.
   */
  protected String checkMainInitOk(String... args) {
    final ByteArrayOutputStream out = new ByteArrayOutputStream();
    final MemoryPrintStream err = new MemoryPrintStream();
    final int rc = mCli.mainInit(args, out, err.printStream());
    assertNotNull(mCli.mMainListener);
    assertEquals("Error: " + err.toString(), "", err.toString());
    assertEquals(0, rc);
    return out.toString();
  }

  /**
   * Runs the main method of the cli class under the assumption that
   * it should run but produce warnings.
   *
   * @param args command line arguments.
   * @return the contents of stderr.
   */
  protected String checkMainInitWarn(String... args) {
    final ByteArrayOutputStream out = new ByteArrayOutputStream();
    final MemoryPrintStream err = new MemoryPrintStream();
    final int rc = mCli.mainInit(args, out, err.printStream());
    assertNotNull(mCli.mMainListener);
    assertTrue(err.toString().length() > 0);
    assertEquals(err.toString(), 0, rc);
    return err.toString();
  }

  /**
   * Runs the main method of the cli class under the assumption that
   * it should fail for some reason.
   *
   * @param args command line arguments.
   * @return the contents of stderr.
   */
  protected String checkMainInitBadFlags(String... args) {
    final MemoryPrintStream err = new MemoryPrintStream();
    assertEquals(1, mCli.mainInit(args, TestUtils.getNullOutputStream(), err.printStream()));
    return err.toString();
  }
}
