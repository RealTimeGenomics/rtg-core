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
package com.rtg.simulation.cnv;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.LogRecord;
import com.rtg.util.io.LogStream;
import com.rtg.util.test.FileHelper;

/**
 * Test the corresponding class
 */
public class CnvSimulatorCliTest extends AbstractCliTest {

  private File mDir;

  @Override
  public void setUp() throws IOException {
    super.setUp();
    mDir = FileHelper.createTempDirectory();
    Diagnostic.setLogStream();
  }

  @Override
  public void tearDown() throws IOException  {
    super.tearDown();
    assertTrue(FileHelper.deleteAll(mDir));
    mDir = null;
  }

  private static final OutputStream NULL_STREAM = TestUtils.getNullOutputStream();
  private static final PrintStream NULL_PRINTSTREAM = TestUtils.getNullPrintStream();

  public void testName() {
    final CnvSimulatorCli sim = new CnvSimulatorCli();
    assertEquals("rtg cnvsim", sim.applicationName() + " " + sim.moduleName());
  }

  /**
   * Test method for {@link CnvSimulatorCli}.
   */
  public final void testGetCFlags() {
    checkHelp(
        "input",
        "output",
        "output file with CNV information",
        "print help on command-line flag usage",
        "approximate minimum percent of template region generated with CNV effects on",
        "number of regions generated with CNV effects on",
        "seed for the random number generator",
    "secondary output SDF");
  }

  public void testFlags() throws IOException {
    final LogStream logStream = new LogRecord();
    Diagnostic.setLogStream(logStream);

    final File out = File.createTempFile("cnvsimulator", "out", mDir);
    final File twin = File.createTempFile("cnvsimulator", "twin", mDir);
    final File cnvsfile = File.createTempFile("cnvsimulator", "test.cnv", mDir);
    assertTrue(out.delete());
    assertTrue(twin.delete());
    assertTrue(cnvsfile.delete());
    final File in = ReaderTestUtils.getDNADir(mDir);
    try {
      assertEquals(1, new CnvSimulatorCli().mainInit(new String[0], NULL_STREAM, NULL_PRINTSTREAM));
      assertEquals(1, new CnvSimulatorCli().mainInit(new String[] {"-i", in.getPath()}, NULL_STREAM, NULL_PRINTSTREAM));
      assertEquals(1, new CnvSimulatorCli().mainInit(new String[] {"-i", in.getPath(), "-o", out.getPath(), "-O", twin.getPath(), "-s", cnvsfile.getPath()}, NULL_STREAM, NULL_PRINTSTREAM));
      } finally {
      assertTrue(FileHelper.deleteAll(in));
      assertTrue(!out.exists() || FileHelper.deleteAll(out));
      assertTrue(!twin.exists() || FileHelper.deleteAll(twin));
      assertTrue(!cnvsfile.exists() || cnvsfile.delete());
    }
  }

  public void testFlagsErrors() throws IOException {
    final File out = File.createTempFile("cnvsimulator", "out", mDir);
    final File twin = File.createTempFile("cnvsimulator", "twin", mDir);
    final File cnvsfile = File.createTempFile("cnvsimulator", "test.cnv", mDir);
    final File in = ReaderTestUtils.getDNADir(mDir);
    try {
      checkFlagsError(new String[] {"--input", new File(mDir, "file is not here").getPath()
              , "--output", out.getPath()
              , "--twin", twin.getPath()
              , "--cnv-file", cnvsfile.getPath()
      }
      , "The specified SDF, \"" + new File(mDir, "file is not here").getPath() + "\", does not exist.");
      checkFlagsError(new String[] {"--input", in.getPath()
              , "--output", out.getPath()
              , "--cnv-file", cnvsfile.getPath()
              , "--twin", twin.getPath()
      }
      , "The directory \"" + out.getPath() + "\" already exists. Please remove it first or choose a different directory.");
      FileHelper.deleteAll(out);
      checkFlagsError(new String[] {"--input", in.getPath()
              , "--output", out.getPath()
              , "--cnv-file", cnvsfile.getPath()
              , "--twin", twin.getPath()
      }
      , "The directory \"" + twin.getPath() + "\" already exists. Please remove it first or choose a different directory.");
      FileHelper.deleteAll(twin);
      checkFlagsError(new String[] {"--input", in.getPath()
              , "--output", out.getPath()
              , "--cnv-file", cnvsfile.getPath()
              , "--twin", twin.getPath()
              , "-Z"
      }
      , "Error: CNV file already exists");
      FileHelper.deleteAll(cnvsfile);
      checkFlagsError(new String[] {"--input", in.getPath()
              , "--output", out.getPath()
              , "--cnv-file", cnvsfile.getPath()
              , "--twin", twin.getPath()
              , "--Xpriors", "monkey_priors_should_be_missing"
      }
      , "Invalid prior option \"monkey_priors_should_be_missing\""/*"Unable to locate \"monkey_priors_should_be_missing\" as a properties file for priors."*/);
      checkFlagsError(new String[] {"--input", in.getPath()
              , "--output", out.getPath()
              , "--cnv-file", cnvsfile.getPath()
              , "--twin", twin.getPath()
              , "--cnv-percent", "110%"
      }
      , "Invalid value \"110%\" for \"--cnv-percent\"");
    } finally {
      //final String str = logStream.toString();
      assertTrue(FileHelper.deleteAll(in));
      assertTrue(!out.exists() || FileHelper.deleteAll(out));
      assertTrue(!twin.exists() || FileHelper.deleteAll(twin));
      assertTrue(!cnvsfile.exists() || cnvsfile.delete());
    }
  }

  public void checkFlagsError(String[] args, String expectedError) {
      final ByteArrayOutputStream baos = new ByteArrayOutputStream();
      final PrintStream err = new PrintStream(baos);
      assertEquals(1, new CnvSimulatorCli().mainInit(args, NULL_STREAM, err));
      err.flush();
      TestUtils.containsAll(baos.toString(), expectedError);
  }

  @Override
  protected AbstractCli getCli() {
    return new CnvSimulatorCli();
  }

}
