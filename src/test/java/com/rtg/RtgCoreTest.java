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
package com.rtg;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;

import com.rtg.launcher.globals.GlobalFlags;
import com.rtg.util.Constants;
import com.rtg.util.TestUtils;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.Talkback;

import junit.framework.TestCase;

/**
 */
public class RtgCoreTest extends TestCase {

  @Override
  public void setUp() {
    Diagnostic.setLogStream();
    GlobalFlags.resetAccessedStatus();
    CommandLine.clearCommandArgs();
  }

  @Override
  public void tearDown() {
    Diagnostic.setLogStream();
    GlobalFlags.resetAccessedStatus();
    CommandLine.clearCommandArgs();
    Talkback.setModuleName(null);
  }

  public void testMessage() {
    final ByteArrayOutputStream bout = new ByteArrayOutputStream();
    final ByteArrayOutputStream berr = new ByteArrayOutputStream();

    try (PrintStream err = new PrintStream(berr)) {
      assertEquals(1, new RtgCore().intMain(new String[]{""}, bout, err));
    } finally {
      try {
        bout.close();
      } catch (final IOException e) {
        // too bad
      }
      //out.close();

      // too bad
      try {
        berr.close();
      } catch (final IOException e) {
      }
    }
    assertEquals("", bout.toString());
    TestUtils.containsAll(berr.toString(), "Usage:", "map", "cgmap", "format", "sdfsplit", "version",
                          "cg2sdf", "sdf2fasta", "sdf2fastq", "snp", "genomesim", "readsim",
                          "sdfstats", "samstats", "sdfsplit", "vcffilter");
  }

  public void testModuleMessage() {
    final ByteArrayOutputStream bout = new ByteArrayOutputStream();
    final ByteArrayOutputStream berr = new ByteArrayOutputStream();

    try (PrintStream err = new PrintStream(berr)) {
      assertEquals(1, new RtgCore().intMain(new String[]{ToolsCommand.FORMAT.getCommandName()}, bout, err));
    } finally {
      try {
        bout.close();
      } catch (final IOException e) {
        // too bad
      }
      //out.close();

      // too bad
      try {
        berr.close();
      } catch (final IOException e) {
      }
    }
    assertEquals("", bout.toString());
    GlobalFlags.resetAccessedStatus();

    final ByteArrayOutputStream busage = new ByteArrayOutputStream();
    final ByteArrayOutputStream out = new ByteArrayOutputStream();
    try (PrintStream usage = new PrintStream(busage)) {
      ToolsCommand.FORMAT.mainInit(new String[0], out, usage);  //spits out expected error from format into usage stream
    } finally {
      try {
        busage.close();
      } catch (final IOException e) {
        // too bad
      }
    }
    assertEquals(busage.toString(), berr.toString());
  }

  private void checkRunningHelp(final String help) throws Exception {
    final ByteArrayOutputStream bout = new ByteArrayOutputStream();
    final ByteArrayOutputStream berr = new ByteArrayOutputStream();

    try (PrintStream err = new PrintStream(berr)) {
      assertEquals(0, new RtgCore().intMain(new String[]{help}, bout, err));
      bout.flush();
    } finally {
      try {
        bout.close();
      } catch (final IOException e) {
        // too bad
      }
      //out.close();

      try {
        berr.close();
      } catch (final IOException e) {
        // too bad
      }
    }
    final ByteArrayOutputStream bout2 = new ByteArrayOutputStream();
    final ByteArrayOutputStream berr2 = new ByteArrayOutputStream();
    try (PrintStream err2 = new PrintStream(berr2)) {
      assertEquals(0, new RtgCore().intMain(new String[]{help, "map"}, bout2, err2));
      err2.flush();
    } finally {
      try {
        bout2.close();
      } catch (final IOException e) {
        // too bad
      }

      try {
        berr2.close();
      } catch (final IOException e) {
        // too bad
      }
    }
    TestUtils.containsAll(bout.toString(), "Usage:", "map", "cgmap", "format", "sdfsplit", "version");
    TestUtils.containsAll(bout2.toString(), "Usage: rtg map", "File Input/Output", "Sensitivity Tuning", "Reporting", "Utility");
    assertEquals("", berr.toString());
    assertEquals("", berr2.toString());
  }

  public void testHelpFlags() throws Exception {
    checkRunningHelp("help");
    checkRunningHelp("-h");
    checkRunningHelp("--help");
    checkRunningHelp("-help");
  }

  private void checkRunningOutputsVersion(final String version, int retCode) throws Exception {
    final ByteArrayOutputStream bout = new ByteArrayOutputStream();

    try (PrintStream err = new PrintStream(bout)) {
      assertEquals(retCode, new RtgCore().intMain(new String[]{version}, bout, err));
      bout.flush();
      assertTrue(bout.toString(), bout.toString().contains("Version: "));
    } finally {
      try {
        bout.close();
      } catch (final IOException e) {
        // too bad
      }
      //out.close();

    }
  }

  public void testVersionFlags() throws Exception {
    checkRunningOutputsVersion("version", 0);
    checkRunningOutputsVersion("--version", 0);
  }

  public void testMessageWithNull() {
    final ByteArrayOutputStream bout = new ByteArrayOutputStream();
    final ByteArrayOutputStream berr = new ByteArrayOutputStream();

    try (PrintStream err = new PrintStream(berr)) {
      assertEquals(1, new RtgCore().intMain(null, bout, err));
    } finally {
      try {
        bout.close();
      } catch (final IOException e) {
        // too bad
      }
      //out.close();

      try {
        berr.close();
      } catch (final IOException e) {
        // too bad
      }
    }

    assertEquals("", bout.toString());
    //assertEquals(USAGE_STR, berr.toString());
  }

  public void testMessageWithInvalidArgument() {
    final ByteArrayOutputStream bout = new ByteArrayOutputStream();
    final ByteArrayOutputStream berr = new ByteArrayOutputStream();

    try (PrintStream err = new PrintStream(berr)) {
      assertEquals(1, new RtgCore().intMain(new String[]{"foo"}, bout, err));
    } finally {
      try {
        bout.close();
      } catch (final IOException e) {
        // too bad
      }
      //out.close();

      try {
        berr.close();
      } catch (final IOException e) {
        // too bad
      }
    }

    TestUtils.containsAll(berr.toString(), "Usage:", "map", "cgmap", "format", "sdfsplit", "version",
        "cg2sdf", "sdf2fasta", "sdf2fastq", "snp", "genomesim", "readsim",
        "sdfstats", "samstats", "sdfsplit", "vcffilter");
  }

  public void testRunningModule() {
    final ByteArrayOutputStream bout = new ByteArrayOutputStream();
    final ByteArrayOutputStream berr = new ByteArrayOutputStream();
    try (PrintStream err = new PrintStream(berr)) {
      assertEquals(1, new RtgCore().intMain(new String[]{"format"}, bout, err));
    } finally {
      try {
        bout.close();
      } catch (final IOException e) {
        // too bad
      }
      //out.close();

      try {
        berr.close();
      } catch (final IOException e) {
        // too bad
      }
    }

    assertEquals("", bout.toString());
    assertTrue(berr.toString().length() > 0);

  }

  public void testShift() {
    final String[] shifted = RtgCore.shift(new String[]{"a", "b"});
    assertEquals(1, shifted.length);
    assertEquals("b", shifted[0]);
  }

  public void testErrorMessage() {
    assertEquals("The " + ToolsCommand.FORMAT.getCommandName() + " command has not been enabled by your current license.\nPlease contact " + Constants.SUPPORT_EMAIL_ADDR + " to have this command licensed.", RtgCore.getErrorMessage(ToolsCommand.FORMAT));
  }
}
