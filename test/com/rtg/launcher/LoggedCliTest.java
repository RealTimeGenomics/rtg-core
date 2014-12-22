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
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;

import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.cli.CFlags;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.SlimException;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.LogStream;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 * Tests corresponding class
 */
public class LoggedCliTest extends TestCase {

  static final String MODULENAME = "loggedCliTest";

  File mDir;

  @Override
  public void setUp() throws IOException {
    mDir = FileUtils.createTempDir("loggedclitest", "blah");
  }

  @Override
  public void tearDown() {
    assertTrue(FileHelper.deleteAll(mDir));
    mDir = null;
  }

  private static final class FakeCli extends LoggedCli {

    private File mDir;
    private int mCode;
    private boolean mThrow;

    public FakeCli(File dir, int code, boolean throwException) {
      mDir = dir;
      mCode = code;
      mThrow = throwException;
    }

    final boolean[] mBools = {false};
    @Override
    protected int mainExec(OutputStream out, LogStream log) {
      mBools[0] = true;
      Diagnostic.warning("blah!");
      if (mThrow) {
        throw new SlimException();
      }
      return mCode;
    }

    @Override
    protected File outputDirectory() {
      return new File(mDir, "new");
    }

    @Override
    protected void initFlags() {
      mFlags = new CFlags();
      mFlags.registerOptional("blah", "nothing really");
    }

    @Override
    public String moduleName() {
      return MODULENAME;
    }
  }

  public void testLogging() throws Exception {
    FakeCli tlcli = new FakeCli(mDir, 0, false);

    tlcli.initFlags();
    assertEquals(0, tlcli.mainExec(TestUtils.getNullOutputStream(), TestUtils.getNullPrintStream()));

    File newFile = new File(mDir, "new");
    assertTrue(newFile.exists());
    assertTrue(tlcli.mBools[0]);

    File logFile = null;
    for (File f : newFile.listFiles()) {
      if (f.toString().equals(newFile.toString() + StringUtils.FS + MODULENAME + ".log")) {
        logFile = f;
      }
    }
    assertNotNull(logFile);
    final String logContents = FileUtils.fileToString(logFile);
    TestUtils.containsAll(logContents
        , "Command line arguments: "
        , "Run Id: "
        , "Finished successfully in "
        , " s."
      );
    assertTrue(new File(newFile, "done").exists());
    final String doneContents = FileUtils.fileToString(new File(newFile, "done"));
    assertTrue(doneContents, doneContents.matches("Finished successfully in \\d+ s\\." + StringUtils.LS));
    assertTrue(new File(newFile, "progress").exists());
    final String progressContents = FileUtils.fileToString(new File(newFile, "progress"));
    TestUtils.containsAll(progressContents
        , "Started"
        , "Finished "
      );
  }

  public void testLoggingFail() throws Exception {
    FakeCli tlcli = new FakeCli(mDir, 1, false);
    tlcli.initFlags();

    assertEquals(1, tlcli.mainExec(TestUtils.getNullOutputStream(), TestUtils.getNullPrintStream()));

    File newFile = new File(mDir, "new");
    assertTrue(newFile.exists());
    assertTrue(tlcli.mBools[0]);

    File logFile = null;
    for (File f : newFile.listFiles()) {
      if (f.toString().equals(newFile.toString() + StringUtils.FS + MODULENAME + ".log")) {
        logFile = f;
      }
    }
    assertNotNull(logFile);
    final String logContents = FileUtils.fileToString(logFile);
    TestUtils.containsAll(logContents
        , "Command line arguments: "
        , "Run failed in "
        , " s."
      );
    assertFalse(new File(newFile, "done").exists());
    assertTrue(new File(newFile, "progress").exists());
    final String progressContents = FileUtils.fileToString(new File(newFile, "progress"));
    TestUtils.containsAll(progressContents
        , "Started"
        , "Run failed"
      );
  }

  public void testLoggingThrown() throws Exception {
    FakeCli tlcli = new FakeCli(mDir, 0, true);
    tlcli.initFlags();

    try {
      tlcli.mainExec(TestUtils.getNullOutputStream(), TestUtils.getNullPrintStream());
      fail();
    } catch (SlimException e) {
      //expected
    }

    File newFile = new File(mDir, "new");
    assertTrue(newFile.exists());
    assertTrue(tlcli.mBools[0]);

    File logFile = null;
    for (File f : newFile.listFiles()) {
      if (f.toString().equals(newFile.toString() + StringUtils.FS + MODULENAME + ".log")) {
        logFile = f;
      }
    }
    assertNotNull(logFile);
    final String logContents = FileUtils.fileToString(logFile);
    TestUtils.containsAll(logContents
        , "Command line arguments: "
        , "Run failed in "
        , " s."
      );
    assertFalse(new File(newFile, "done").exists());
    assertTrue(new File(newFile, "progress").exists());
    final String progressContents = FileUtils.fileToString(new File(newFile, "progress"));
    TestUtils.containsAll(progressContents
        , "Started"
        , "Run failed"
      );
  }

  public void testCleanDir() throws Exception {
    Diagnostic.setLogStream();
    final String moduleName = "loggedCliTest";

    final File dir = FileUtils.createTempDir("loggedclitest", "blah");

    try {
      LoggedCli tlcli = new LoggedCli() {

        @Override
        protected int mainExec(OutputStream out, LogStream log) {
          assertTrue(outputDirectory().exists());
          cleanDirectory(); // Signal that we want to clean the output directory (if we created it)
          return 0;
        }

        @Override
        protected File outputDirectory() {
          return new File(dir, "new");
        }

        @Override
        protected void initFlags() {
          mFlags = new CFlags();
        }

        @Override
        public String moduleName() {
          return moduleName;
        }

      };

      // Can we create a directory
      ByteArrayOutputStream errbaos = new ByteArrayOutputStream();
      tlcli.createDirectory(tlcli.outputDirectory());
      tlcli.cleanDirectory();
      tlcli.initFlags();
      File newFile = new File(dir, "new");
      assertTrue(newFile.exists());

      // Output directory now already exists. We should not delete it when we run
      assertEquals(0, tlcli.mainExec(new ByteArrayOutputStream(), new PrintStream(errbaos)));
      assertTrue(newFile.exists());

      // If we use an output directory that we create, we should delete it
      assertTrue(FileHelper.deleteAll(newFile));
      assertEquals(0, tlcli.mainExec(new ByteArrayOutputStream(), new PrintStream(errbaos)));
      assertFalse(newFile.exists());
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  public void testTimeDifference() {
    assertEquals(1, LoggedCli.timeDifference(2000, 1000), 0);
  }
}
