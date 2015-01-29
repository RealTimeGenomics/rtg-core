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

import static com.rtg.util.StringUtils.LS;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;

import com.rtg.util.IORunnable;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.Params;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.CliDiagnosticListener;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.LogRecord;
import com.rtg.util.io.LogStream;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 */
public class ParamsCliTest extends TestCase {

  protected File mDir;

  @Override
  public void setUp() throws IOException {
    Diagnostic.setLogStream();
    GlobalFlags.resetAccessedStatus();
    mDir = FileHelper.createTempDirectory();
  }

  @Override
  public void tearDown() {
    assertTrue(FileHelper.deleteAll(mDir));
    mDir = null;
  }

  public void test() throws IOException, InvalidParamsException {
    final File outDir = FileHelper.createTempDirectory();
    final ParamsCli<MockCliParams> cli = new MockCli(outDir);
    try {
      final ByteArrayOutputStream out = new ByteArrayOutputStream();
      assertEquals(1, cli.mainInit(new String[0], out, StreamUtil.make()));
      assertEquals("Mock task did something", new String(out.toByteArray()));
    } finally {
      assertTrue(FileHelper.deleteAll(outDir));
    }
  }

  public void testParamsErr() throws IOException, InvalidParamsException {
    final File subject = BuildTestUtils.prereadDNA(mDir, ">s" + LS + "ACTGA" + LS);
    final String su = subject.getAbsolutePath();
    final File query = BuildTestUtils.prereadDNA(mDir, ">q" + LS + "CTGA" + LS);
    final String qu = query.getAbsolutePath();
    final File dir = FileUtils.createTempDir("test", "outdir");
    FileHelper.deleteAll(dir);
    final String di = dir.getAbsolutePath();
    final ByteArrayOutputStream out = new ByteArrayOutputStream();
    final PrintStream err = StreamUtil.make();
    final ParamsCli<MockCliParams> cli = new MockCli();
    assertEquals(1, cli.mainInit(new String[] {"-o", di, "-i", su, "-x", qu, "-l", "0"}, out, err));

    assertFalse(dir.exists());
    assertEquals("", new String(out.toByteArray()));
    final String errExp = ""
      + "Error: Unknown flag -o" + LS + LS
      + "Usage: MockCli [OPTION]..." + LS
      + LS
      + "Try '--help' for more information" + LS
      ;
    final String foo = err.toString();
    assertEquals(errExp, foo);
    //System.err.println(dlog);

    assertTrue(FileHelper.deleteAll(query));
    assertTrue(FileHelper.deleteAll(subject));
  }

  public final void testSomethingOK() throws IOException, InvalidParamsException {
    final ParamsCli<?> cli = new MockCli();
    final ByteArrayOutputStream berr = new ByteArrayOutputStream();
    final PrintStream err = new PrintStream(berr);
    final ByteArrayOutputStream bout = new ByteArrayOutputStream();
    final PrintStream out = new PrintStream(bout);
    cli.createFlags(out, err);
    cli.initFlags();
    cli.mFlags.setFlags();
    final Params p = cli.makeParams();
    assertNotNull(p);
    final String estr = berr.toString();
    //System.err.println(estr);
    assertEquals("", estr);
  }

  //error in validator
  public final void testErrV() throws IOException, InvalidParamsException {
    final LogStream log = new LogRecord();
    Diagnostic.setLogStream(log);
    final ByteArrayOutputStream berr = new ByteArrayOutputStream();
    final PrintStream err = new PrintStream(berr);
    final ByteArrayOutputStream bout = new ByteArrayOutputStream();
    final PrintStream out = new PrintStream(bout);
    final ParamsCli<?> cli = new MockCli();
    cli.createFlags(out, err);
    cli.initFlags();
    try {
      final CliDiagnosticListener listener = new CliDiagnosticListener(err);
      Diagnostic.addListener(listener);
      try {
        assertFalse(cli.handleFlags(new String[] {"-v"}, out, err));
        //System.err.println("err:" + estr);
        //System.err.println("log:" + logs);
      } finally {
        Diagnostic.removeListener(listener);
      }
    } finally {
      Diagnostic.setLogStream();
      err.close();
    }
    final String exp = "The specified flag \"maxgap\" has invalid value \"42\". It should be greater than or equal to 1.";
    final String estr = berr.toString();
    final String logs = log.toString();
    assertTrue(logs.contains(exp));
    assertTrue(estr.contains(exp));
  }

  //error in constructor
  public final void testErrC() throws IOException, InvalidParamsException {
    final LogStream log = new LogRecord();
    Diagnostic.setLogStream(log);
    final ByteArrayOutputStream berr = new ByteArrayOutputStream();
    final PrintStream err = new PrintStream(berr);
    final ByteArrayOutputStream bout = new ByteArrayOutputStream();
    final PrintStream out = new PrintStream(bout);
    final ParamsCli<?> cli = new MockCli();
    cli.createFlags(out, err);
    cli.initFlags();
    try {
      final CliDiagnosticListener listener = new CliDiagnosticListener(err);
      Diagnostic.addListener(listener);
      try {
        cli.mFlags.setFlags("-c");
        try {
          cli.makeParams();
          fail();
        } catch (final InvalidParamsException e) {
          assertEquals(ErrorType.INVALID_MAX_INTEGER_FLAG_VALUE, e.getErrorType());
        }
        //System.err.println("err:" + estr);
        //System.err.println("log:" + logs);
      } finally {
        Diagnostic.removeListener(listener);
      }
    } finally {
      Diagnostic.setLogStream();
      err.close();
    }
  }

  public void testGood() throws IOException {
    final File dir = FileUtils.createTempDir("paramscli", "test");
    try {
      final MemoryPrintStream err = new MemoryPrintStream();
      TestParamsCli cli = new TestParamsCli(false, false, dir);
      final int code = cli.mainInit(new String[] {}, TestUtils.getNullOutputStream(), err.printStream());
      assertEquals(err.toString(), 0, code);
      assertTrue(cli.hasRan());
      assertTrue(dir.exists());
      assertTrue(FileHelper.deleteAll(dir));
      assertEquals(0, cli.mainInit(new String[] {}, TestUtils.getNullOutputStream(), err.printStream()));
      assertTrue(dir.exists());
    } finally {
      assertTrue(!dir.exists() || FileHelper.deleteAll(dir));
    }
  }

  public void testMakeParamsHandling() throws IOException {
    final File dir = FileUtils.createTempDir("paramscli", "test");
    try {
      final MemoryPrintStream err = new MemoryPrintStream();
      TestParamsCli cli = new TestParamsCli(true, false, dir);
      final int code = cli.mainInit(new String[] {}, TestUtils.getNullOutputStream(), err.printStream());
      assertEquals(err.toString(), 1, code);
      assertTrue(err.toString().contains("Test error"));
      assertFalse(cli.hasRan());
      assertTrue(dir.exists());
      assertTrue(FileHelper.deleteAll(dir));
      assertEquals(1, cli.mainInit(new String[] {}, TestUtils.getNullOutputStream(), err.printStream()));
      assertFalse(dir.exists());
    } finally {
      assertTrue(!dir.exists() || FileHelper.deleteAll(dir));
    }
  }

  public void testMakeParamsHandling2() throws IOException {
    final File dir = FileUtils.createTempDir("paramscli", "test");
    try {
      final MemoryPrintStream err = new MemoryPrintStream();
      TestParamsCli cli = new TestParamsCli(true, true, dir);
      final int code = cli.mainInit(new String[] {}, TestUtils.getNullOutputStream(), err.printStream());
      assertEquals(err.toString(), 1, code);
      assertTrue(err.toString().contains("Test error"));
      assertFalse(cli.hasRan());
      assertTrue(dir.exists());
      assertTrue(FileHelper.deleteAll(dir));
      assertEquals(1, cli.mainInit(new String[] {}, TestUtils.getNullOutputStream(), err.printStream()));
      assertFalse(dir.exists());
    } finally {
      assertTrue(!dir.exists() || FileHelper.deleteAll(dir));
    }
  }

  public void testCheckFiles() throws IOException, InvalidParamsException {
    final File dir = FileUtils.createTempDir("paramscli", "test");
    try {
      ArrayList<Object> files = new ArrayList<>();
      final File created = new File(dir, "exists");
      final File fobar = new File(dir, "fobar");
      assertTrue(created.createNewFile());
      files.add(fobar);
      files.add(created);
      files.add(created);
      final MemoryPrintStream err = new MemoryPrintStream();
      Diagnostic.setLogStream(err.printStream());
      try {
        ParamsCli.checkFiles(files);
        fail();
      } catch (final InvalidParamsException e) {
        assertEquals(ErrorType.INFO_ERROR, e.getErrorType());
      } finally {
        Diagnostic.setLogStream();
      }
      assertTrue(err.toString().contains("File not found"));
//      assertTrue(err.toStringNoFlush().contains("1 specified files were not found"));
      files.remove(fobar);
      final Collection<File> collfiles = ParamsCli.checkFiles(files);
      assertEquals(1, collfiles.size());
      assertTrue(collfiles.contains(created));
    } finally {
      assertTrue(!dir.exists() || FileHelper.deleteAll(dir));
    }
  }

  private static class BogusParams implements Params {
    @Override
    public void close() {
    }
  }
  private static class TestParamsCli extends ParamsCli<BogusParams> {
    private final boolean mThrowException;
    private final boolean mNoTalkback;
    private boolean mRan;
    private final File mDir;

    public TestParamsCli(boolean throwException, boolean notalkback, File dir) {
      mThrowException = throwException;
      mNoTalkback = notalkback;
      mDir = dir;
    }

    public boolean hasRan() {
      return mRan;
    }
    @Override
    protected BogusParams makeParams() throws InvalidParamsException, IOException {
      if (mThrowException) {
        if (mNoTalkback) {
          throw new NoTalkbackSlimException(ErrorType.INFO_ERROR, "Test error");
        }
        throw new InvalidParamsException(ErrorType.INFO_ERROR, "Test error");
      } else {
        return new BogusParams();
      }
    }
    @Override
    protected IORunnable task(BogusParams params, OutputStream out) {
      return new IORunnable() {

        @Override
        public void run() {
          mRan = true;
        }
      };
    }
    @Override
    protected File outputDirectory() {
      return mDir;
    }
    @Override
    protected void initFlags() {
    }
    @Override
    public String moduleName() {
      return "test";
    }
  }
}
