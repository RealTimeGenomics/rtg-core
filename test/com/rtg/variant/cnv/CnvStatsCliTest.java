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
package com.rtg.variant.cnv;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.launcher.BuildTestUtils;
import com.rtg.simulation.cnv.CnvSimulatorCli;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.TestCFlags;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.LogRecord;
import com.rtg.util.io.LogStream;
import com.rtg.util.test.FileHelper;

/**
 * Test the corresponding class
 */
public class CnvStatsCliTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new CnvStatsCli();
  }

  static final String LS = StringUtils.LS;
  static final String TB = "\t";
  private File mDir;

  @Override
  public void setUp() throws IOException {
    super.setUp();
    mDir = FileHelper.createTempDirectory();
    Diagnostic.setLogStream();
  }

  @Override
  public void tearDown() throws IOException {
    super.tearDown();
    assertTrue(FileHelper.deleteAll(mDir));
    mDir = null;
  }

  private static final OutputStream NULL_STREAM = TestUtils.getNullOutputStream();
  private static final PrintStream NULL_PRINTSTREAM = TestUtils.getNullPrintStream();


//  public void testName() {
//    final CnvStatsCli sim = new CnvStatsCli();
//    assertEquals("rtg cnvstats", sim.applicationName() + " " + sim.moduleName());
//  }

  //private static final String USAGE = " rtg cnvsim [OPTION]... -n INT -s FILE -i SDF -o SDF -O SDF";


  /**
   * Test method for {@link CnvSimulatorCli}.
   */
  public final void testGetCFlags() {
    final CFlags flags = new CFlags();
    final CnvStatsCli stats = new CnvStatsCli();
    stats.initFlags(flags);

    TestCFlags.check(flags,
        "CNV generation file",
        "CNV detection file",
    "output directory for results");
  }

  public void checkErrorMessage(final String[] args0, final String exp) throws InvalidParamsException, IOException {
    final LogStream logStream = new LogRecord();
    Diagnostic.setLogStream(logStream);
    final ByteArrayOutputStream berr = new ByteArrayOutputStream();
    final PrintStream err = new PrintStream(berr);
    final File subjectsDir = FileHelper.createTempDirectory();
    try {
      final CnvStatsCli cnvstats = (CnvStatsCli) mCli;
        checkHandleFlagsErr(args0);
        try {
          assertEquals(1, cnvstats.mainInit(args0, NULL_STREAM,  err));
          err.flush();
          //System.err.println("err: " + berr.toString());
          //err.flush();
          TestUtils.containsAll(berr.toString(), exp);
        } catch (final Exception e) {
          //assertEquals(null, e.getMessage()); // Cannot assume this on .NET
          final String str = logStream.toString();
          //System.err.println(str);
          assertTrue("Exception:" + e.getClass().getName() + " Actual: " + str + "\n" + "Expected to contain: " + exp, str.contains(exp));
        }
    } finally {
      assertTrue(FileHelper.deleteAll(subjectsDir));
    }
  }

  public void testFlags() throws IOException {
    final LogStream logStream = new LogRecord();
    Diagnostic.setLogStream(logStream);

    final File cnvsfile1 = File.createTempFile("cnvsimulator", "test1.cnv", mDir);
    final File cnvsfile2 = File.createTempFile("cnvsimulator", "test2.cnv", mDir);
    assertTrue(cnvsfile1.delete());
    assertTrue(cnvsfile2.delete());
    final File in = BuildTestUtils.prereadDNA(mDir, ">a\nacgtacgatcagcatctgac\n");
    try {
      assertEquals(1, new CnvSimulatorCli().mainInit(new String[0], NULL_STREAM, NULL_PRINTSTREAM));
      assertEquals(1, new CnvSimulatorCli().mainInit(new String[] {"-g", cnvsfile1.getPath()}, NULL_STREAM, NULL_PRINTSTREAM));
      assertEquals(1, new CnvSimulatorCli().mainInit(new String[] {"-s", cnvsfile2.getPath()}, NULL_STREAM, NULL_PRINTSTREAM));
    } finally {
      Diagnostic.setLogStream();
      //final String str = logStream.toString();
      assertTrue(FileHelper.deleteAll(in));
      assertTrue(!cnvsfile1.exists() || FileHelper.deleteAll(cnvsfile1));
      assertTrue(!cnvsfile2.exists() || FileHelper.deleteAll(cnvsfile2));
    }
  }
  static final String SIMULATED_1 = "#Version v2.0blabla " + LS
    + "#CL     cnvsim blabla " + LS
    + " #Seq" + TB + "start" + TB + "end" + TB + "label" + TB + "cn" + TB + "bp-cn" + TB + "error" + LS
    + " seq1" + TB + "0" + TB + "39992" + TB + "cnv" + TB + "2" + TB + "0" + TB + "0.0" + LS;

  public void testFlagErrorMessages() throws IOException, InvalidParamsException {
    final LogStream logStream = new LogRecord();
    Diagnostic.setLogStream(logStream);

    final File cnvsfileSimulated = File.createTempFile("cnvsimulator", "test1.cnv", mDir);
    FileUtils.stringToFile(SIMULATED_1, cnvsfileSimulated);
    final File cnvsfileGenerated = File.createTempFile("cnvsimulator", "test1.cnv", mDir);
    FileUtils.stringToFile(SIMULATED_1, cnvsfileGenerated);

    final File doesnotexist = File.createTempFile("cnvsimulator", "test1.cnv", mDir);
    assertTrue(doesnotexist.delete());
    final File in = BuildTestUtils.prereadDNA(mDir, ">a\nacgtacgatcagcatctgac\n");
    try {
      checkErrorMessage(new String[] {"-s", cnvsfileSimulated.getPath(), "-g"}, "Expecting value for flag -g");
      checkErrorMessage(new String[] {"-s", cnvsfileSimulated.getPath(), "-g", doesnotexist.getPath()},
      "CNV generation file doesn't exist");
      checkErrorMessage(new String[] {"-s", doesnotexist.getPath(), "-g", cnvsfileGenerated.getPath()},
      "CNV detection file doesn't exist");
    } finally {
      Diagnostic.setLogStream();
      //final String str = logStream.toString();
      assertTrue(FileHelper.deleteAll(in));
      assertTrue(!cnvsfileSimulated.exists() || FileHelper.deleteAll(cnvsfileSimulated));
      assertTrue(!cnvsfileGenerated.exists() || FileHelper.deleteAll(cnvsfileGenerated));
    }
  }
}
