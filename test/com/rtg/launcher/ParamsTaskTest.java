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
import java.io.StringWriter;

import com.rtg.util.InvalidParamsException;
import com.rtg.util.cli.CFlags;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.LogRecord;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 */
public class ParamsTaskTest extends TestCase {

  protected File mDir;

  @Override
  public void setUp() throws IOException {
    Diagnostic.setLogStream();
    mDir = FileHelper.createTempDirectory();
  }

  @Override
  public void tearDown() {
    assertTrue(FileHelper.deleteAll(mDir));
    mDir = null;
  }

  public void test() throws IOException, InvalidParamsException {
    final LogRecord log = new LogRecord();
    Diagnostic.setLogStream(log);
    try {
      final File subject = BuildTestUtils.prereadDNA(mDir, ">s" + LS + "ACTGA" + LS);
      final String su = subject.getAbsolutePath();
      final File query = BuildTestUtils.prereadDNA(mDir, ">q" + LS + "CTGA" + LS);
      final String qu = query.getAbsolutePath();
      final File dir = FileUtils.createTempDir("test", "outdir");
      final String di = dir.getAbsolutePath();
      final MockCliParams pr = getParamsStatic(new String[]{"-o", di, "-i", su, "-x", qu, "-l", "4", "-Z"});
      final File diro = FileHelper.createTempDirectory();
      pr.setDirectory(diro);
      final ByteArrayOutputStream out = new ByteArrayOutputStream();
      final ParamsTask<MockCliParams, NoStatistics> pc = new DummyTask(pr, out);
      assertEquals(pr, pc.parameters());

      pc.run();

      assertEquals("", out.toString());
      assertTrue(FileHelper.deleteAll(diro));
      assertTrue(FileHelper.deleteAll(query));
      assertTrue(FileHelper.deleteAll(subject));
      assertTrue(FileHelper.deleteAll(dir));
    } finally {
      Diagnostic.setLogStream();
    }
  }

  private class DummyTask extends MockTask {
    DummyTask(final MockCliParams params, final OutputStream defaultOutput) {
      super(params, defaultOutput);
    }
  }

  static MockCliParams getParamsStatic(final String[] args) throws InvalidParamsException, IOException {
    final Appendable out = new StringWriter();
    final Appendable err = new StringWriter();
    final CFlags flags = new CFlags("testMockParams", out, err);
    //MockCliParams.initFlags(flags);
    flags.setFlags(args);
    return new MockCliParams(flags);
  }

}
