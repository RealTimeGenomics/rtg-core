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

import java.io.File;
import java.io.IOException;

import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.MemoryPrintStream;

import junit.framework.TestCase;

/**
 */
public class DummyStatisticsTest extends TestCase {

  private static class DummyStatistics extends AbstractStatistics {

    protected String mStats = null;

    /**
     * @param outputDirectory The base output directory to generate statistics and reports in. May be null if no statistics or reports are to be generated.
     */
    public DummyStatistics(File outputDirectory) {
      super(outputDirectory);
    }

    @Override
    protected String getStatistics() {
      return mStats;
    }

    @Override
    public void generateReport() { }
  }

  @Override
  public void tearDown() {
    Diagnostic.setLogStream();
  }

  public void test() throws IOException {
    final MemoryPrintStream ps = new MemoryPrintStream();
    final MemoryPrintStream ps2 = new MemoryPrintStream();

    Diagnostic.setLogStream(ps2.printStream());
    final DummyStatistics stat = new DummyStatistics(null);
    stat.printStatistics(null);
    assertEquals("", ps.toString());
    stat.printStatistics(null);
    assertEquals("", ps.toString());
    stat.mStats = "Line1" + StringUtils.LS + "Line2";

    stat.printStatistics(ps.outputStream());
    assertEquals("Line1" + StringUtils.LS + "Line2", ps.toString());
    assertEquals(2, ps.toString().split(StringUtils.LS).length);

    TestUtils.containsAll(ps2.toString(), "Line1", "Line2"); //logged statistics
    assertEquals(2, ps2.toString().split(StringUtils.LS).length);
  }
}
