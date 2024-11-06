/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
      "secondary output SDF"
    );
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
      assertEquals(0, new CnvSimulatorCli().mainInit(new String[] {"-i", in.getPath(), "-o", out.getPath(), "-O", twin.getPath(), "-s", cnvsfile.getPath()}, NULL_STREAM, NULL_PRINTSTREAM));
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
      , "Invalid value \"110%\" for flag --cnv-percent");
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
