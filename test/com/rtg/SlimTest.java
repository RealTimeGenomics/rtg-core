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
package com.rtg;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;

import com.rtg.launcher.GlobalFlags;
import com.rtg.util.Constants;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.Talkback;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * Class to test <code>Slim</code>
 *
 */
public class SlimTest extends TestCase {

  public SlimTest(final String name) {
    super(name);
  }

  public static Test suite() {
    return new TestSuite(SlimTest.class);
  }
  /**
   * Main to run from tests from command line.
   * @param args ignored.
   */
  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }

  @Override
  public void setUp() {
    Diagnostic.setLogStream();
    GlobalFlags.resetAccessedStatus();
  }

  public void testMessage() {
    final ByteArrayOutputStream bout = new ByteArrayOutputStream();
    final ByteArrayOutputStream berr = new ByteArrayOutputStream();

    try (PrintStream err = new PrintStream(berr)) {
      assertEquals(1, new Slim().intMain(new String[]{""}, bout, err));
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
      assertEquals(1, new Slim().intMain(new String[]{ToolsCommand.FORMAT.getCommandName()}, bout, err));
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
      assertEquals(0, new Slim().intMain(new String[]{help}, bout, err));
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
      assertEquals(0, new Slim().intMain(new String[]{help, "map"}, bout2, err2));
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
      assertEquals(retCode, new Slim().intMain(new String[]{version}, bout, err));
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
      assertEquals(1, new Slim().intMain(null, bout, err));
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
      assertEquals(1, new Slim().intMain(new String[]{"foo"}, bout, err));
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
      assertEquals(1, new Slim().intMain(new String[]{"format"}, bout, err));
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

  public void testSearch() {
    final ByteArrayOutputStream bout = new ByteArrayOutputStream();
    final ByteArrayOutputStream berr = new ByteArrayOutputStream();
    try (PrintStream err = new PrintStream(berr)) {
      assertEquals(1, new Slim().intMain(new String[]{"format"}, bout, err));
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
    final String[] shifted = Slim.shift(new String[]{"a", "b"});
    assertEquals(1, shifted.length);
    assertEquals("b", shifted[0]);
  }

  public void testFindXlog() {
    final PrintStream systemErr = System.err;
    Object[] objs = Slim.findXlog(new String[]{"--", "--Xlog"});
    assertEquals(2, objs.length);
    assertNull(objs[0]);
    assertTrue(objs[1] instanceof String[]);
    assertEquals(2, ((String[]) objs[1]).length);
    assertEquals("--", ((String[]) objs[1])[0]);
    assertEquals("--Xlog", ((String[]) objs[1])[1]);
    objs = Slim.findXlog(new String[]{"abc", "--Xlog=blah", "123"});
    assertEquals(2, objs.length);
    assertTrue(objs[0] instanceof String);
    assertEquals("blah", (String) objs[0]);
    assertTrue(objs[1] instanceof String[]);
    assertEquals(2, ((String[]) objs[1]).length);
    assertEquals("abc", ((String[]) objs[1])[0]);
    assertEquals("123", ((String[]) objs[1])[1]);
    objs = Slim.findXlog(new String[]{"123", "--Xlog", "blarg", "abc"});
    assertEquals(2, objs.length);
    assertTrue(objs[0] instanceof String);
    assertEquals("blarg", (String) objs[0]);
    assertTrue(objs[1] instanceof String[]);
    assertEquals(2, ((String[]) objs[1]).length);
    assertEquals("123", ((String[]) objs[1])[0]);
    assertEquals("abc", ((String[]) objs[1])[1]);
    final ByteArrayOutputStream bos = new ByteArrayOutputStream();
    final PrintStream err = new PrintStream(bos);
    System.setErr(err);
    try {
      Slim.findXlog(new String[]{"a", "--Xlog"});
      fail();
    } catch (RuntimeException e) { }
    err.flush();
    TestUtils.containsAll(bos.toString(), "Expected URL after --Xlog");
    System.setErr(systemErr);
    err.close();
  }

  public void testErrorMessage() {
    assertEquals("The " + ToolsCommand.FORMAT.getCommandName() + " command has not been enabled by your current license.\nPlease contact " + Constants.SUPPORT_EMAIL_ADDR + " to have this command licensed.", Slim.getErrorMessage(ToolsCommand.FORMAT));
  }

  @Override
  public void tearDown() {
    // clear the module name so later tests don't report SlimException to the Talkback system
    Talkback.setModuleName(null);
  }
}
