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
package com.rtg.util;

import static com.rtg.util.StringUtils.FS;

import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.io.PrintStream;
import java.io.StringWriter;

import com.rtg.launcher.BuildCommon;
import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.OutputParams;
import com.rtg.launcher.OutputParamsTest;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.TestCFlags;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 */
public class ObjectParamsTest extends TestCase {

  public static Test suite() {
    final TestSuite suite = new TestSuite();
    suite.addTestSuite(ObjectParamsTest.class);
    suite.addTestSuite(OutputParamsTest.class);
    return suite;
  }

  private static final class MockParams extends OutputParams {

    private final Integer mMock;

    MockParams(final CFlags flags, final int mock) {
      super((File) flags.getValue(CommonFlags.OUTPUT_FLAG), flags.isSet(BuildCommon.PROGRESS_FLAG), !flags.isSet(CommonFlags.NO_GZIP));
      mMock = mock;
      append(new Object[] {mMock});
    }

  }

  /**
   * Create a parameter object to be tested.
   * @param args used for construction.
   * @return parameters object.
   */
  protected OutputParams getParams(final String[] args) {
    final CFlags flags = new CFlags("testOutputParams", new StringWriter(), null);
    OutputParams.initFlags(flags);
    flags.setFlags(args);
    final OutputParams params = new MockParams(flags, 1);
    assertTrue(params.closed());
    assertNotNull(params.objects());
    return params;
  }

  public void test() throws Exception {
    final File dirName = FileUtils.createTempDir("output", "");
    assertTrue(dirName.delete());
    final OutputParams a1 = getParams(new String[] {"-o", dirName.getPath(), "-Z"});
    final OutputParams a2 = getParams(new String[] {"-o", dirName.getPath(), "-Z"});
    final OutputParams c = getParams(new String[] {"-o", dirName.getPath()});
    final OutputParams d = getParams(new String[] {"-o", dirName.getPath(), "-P", "-Z"});
    TestUtils.equalsHashTest(new OutputParams[][] {{a1, a2}, {c}, {d}});

    assertEquals("OutputParams output directory=" + dirName.getPath() + " progress=" + Boolean.FALSE.toString() + " zip=" + Boolean.FALSE.toString(), a1.toString());
    assertEquals(dirName.getPath(), a1.directory().toString());
    assertTrue(a1.file("child").toString().endsWith(dirName.getPath() + FS + "child"));
    assertEquals(false, a1.progress());

    assertEquals("OutputParams output directory=" + dirName.getPath() + " progress=" + Boolean.FALSE.toString() + " zip=" + Boolean.TRUE.toString(), c.toString());
    assertEquals(dirName.getPath(), c.directory().toString());
    assertTrue(c.file("child").toString().endsWith(dirName.getPath() + FS + "child"));
    assertEquals(false, c.progress());

    assertEquals("OutputParams output directory=" + dirName.getPath() + " progress=" + Boolean.TRUE.toString() + " zip=" + Boolean.FALSE.toString(), d.toString());
    assertEquals(dirName.getPath(), d.directory().toString());
    assertTrue(d.file("child").toString().endsWith(dirName.getPath() + FS + "child"));
    assertEquals(true, d.progress());

    a1.close();
    a2.close();
    c.close();
    d.close();
  }

  public void testStreamO() throws Exception {
    final File dirName = FileUtils.createTempDir("output", "");
    assertTrue(dirName.delete());
    final OutputParams a = getParams(new String[] {"-o", dirName.getPath(), "-Z"});
    try {
      final File dir = a.directory();
      final File file = new File(dir, "out");
      final PrintStream ps = new PrintStream(a.outStream("out"));
      ps.append("foobar");
      ps.close();
      try (InputStream in = new FileInputStream(file)) {
        assertEquals("foobar", FileUtils.streamToString(in));
      }
    } finally {
      a.close();
      FileHelper.deleteAll(a.directory());
    }
  }

  public void testStreamOZ() throws Exception {
    final File dirName = FileUtils.createTempDir("output", "");
    assertTrue(dirName.delete());
    final OutputParams a = getParams(new String[] {"-o", dirName.getPath()});
    try {
      final File dir = a.directory();
      final File file = new File(dir, "out.gz");
      final PrintStream ps = new PrintStream(a.outStream("out"));
      ps.append("foobar");
      ps.close();
      a.close();

      assertEquals("foobar", FileHelper.zipFileToString(file));
    } finally {
      a.close();
      FileHelper.deleteAll(a.directory());
    }
  }

  public final void testFlags() {
    final CFlags flags = new CFlags("testOutputParams", new StringWriter(), null);
    OutputParams.initFlags(flags);
    TestCFlags.check(flags, "[OPTION]... -o DIR",
                     "directory for output",
                     "do not gzip the output",
                     "report progress"
                     );
  }
}
