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
package com.rtg.graph;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

import com.reeltwo.plot.Graph2D;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;

import junit.framework.TestCase;

/**
 */
public class ParseRocFileTest extends TestCase {

  private static final String ROC = ""
          + "#total baseline variants: 3092754" + StringUtils.LS
          + "#score\ttrue_positives\tfalse_positives" + StringUtils.LS
          + "3.300\t0.000\t15" + StringUtils.LS
          + "2.261\t70000.000\t137" + StringUtils.LS
          + "1.226\t180000.000\t516" + StringUtils.LS
          + "0.700\t406000.000\t11337" + StringUtils.LS
          + "0.533\t1971000.000\t1446920" + StringUtils.LS
          + "0.333\t2071000.000\t1646920" + StringUtils.LS
          + "0.200\t2995295.000\t1864591" + StringUtils.LS;

  public void test() throws IOException {
    try (final TestDirectory dir = new TestDirectory()) {
      final File roc = FileUtils.stringToFile(ROC, new File(dir, "roc.tsv"));
      try (final BufferedInputStream buf = new BufferedInputStream(new FileInputStream(roc))) {
        final DataBundle bundle = ParseRocFile.loadStream(new MyProgressDelegate(), buf, "Monkey", false);

        assertEquals("Monkey", bundle.getTitle());
        assertEquals(3092754, bundle.getTotalVariants());
        assertEquals(1864591.0, bundle.getPlot(1, 7342).getHi(Graph2D.X), 1e-9);
        assertEquals(0.0, bundle.getPlot(1, 7342).getLo(Graph2D.X), 1e-9);
        assertEquals(2995295.0, bundle.getPlot(1, 7342).getHi(Graph2D.Y), 1e-9);
        assertEquals(0.0, bundle.getPlot(1, 7342).getLo(Graph2D.Y), 1e-9);
        assertEquals(8, bundle.getPlot(1, 1).getData().length);
      }
    }
  }

  private static final String ROC_MAL1 = ""
          + "3.300\t0.000" + StringUtils.LS;
  private static final String ROC_MAL2 = ""
          + "3.300\tROBOT\t15" + StringUtils.LS;
  private static final String ROC_MAL3 = ""
          + "3.300\t0.000\tROBOT" + StringUtils.LS;

  private static final String[] MALFORMED = {ROC_MAL1, ROC_MAL2, ROC_MAL3};
  public void testMalformed() throws IOException {
    Diagnostic.setLogStream();
    try (final TestDirectory dir = new TestDirectory()) {
      for (String mal : MALFORMED) {
        final File roc = FileUtils.stringToFile(mal, new File(dir, "roc.tsv"));
        try (final BufferedInputStream buf = new BufferedInputStream(new FileInputStream(roc))) {
          try {
            ParseRocFile.loadStream(new MyProgressDelegate(), buf, "Monkey", false);
            fail("Should get exception: " + mal);
          } catch (RuntimeException | IOException e) {

          }
        }
        assertTrue(roc.delete());
      }
    }
  }

  public void testProgress() throws IOException {
    Diagnostic.setLogStream();
    try (final TestDirectory dir = new TestDirectory()) {
      final StringBuilder sb = new StringBuilder();
      for (int i = 0; i < 500; i++) {
        sb.append("1.0\t5.0\t2").append(StringUtils.LS);
      }
      for (int i = 0; i < 500; i++) {
        sb.append("1.0\t7.0\t9").append(StringUtils.LS);
      }

      final File roc = FileUtils.stringToFile(sb.toString(), new File(dir, "roc.tsv"));
      try (final BufferedInputStream buf = new BufferedInputStream(new FileInputStream(roc))) {
        final MyProgressDelegate progressBarDelegate = new MyProgressDelegate();
        final DataBundle bundle = ParseRocFile.loadStream(progressBarDelegate, buf, "Monkey", true);
        assertEquals(1000, progressBarDelegate.mNumberLines);
        assertEquals(Arrays.asList(new Integer[] {100, 200, 300, 400, 500, 600, 700, 800, 900, 1000}), progressBarDelegate.mProgressLog);
        assertEquals(2, bundle.getPlot(1, 1).getData().length);
        assertEquals(-1, bundle.getTotalVariants());
      }
    }
  }

  public static void main(String[] args) throws Exception {
    try (TestDirectory dir = new TestDirectory()) {
      final File roc = FileUtils.stringToFile(ROC_MAL1, new File(dir, "roc.tsv"));
      final MemoryPrintStream mps = new MemoryPrintStream();
      new RocPlotCli().mainInit(new String[] {roc.getPath()}, mps.outputStream(), mps.printStream());
    }
  }

  private static class MyProgressDelegate implements ProgressDelegate {
    final ArrayList<Integer> mProgressLog = new ArrayList<>();
    int mNumberLines = 0;

    @Override
    public void setProgress(int progress) {
      mProgressLog.add(progress);
    }

    @Override
    public void addFile(int numberLines) {
      mNumberLines += numberLines;
    }

    @Override
    public void done() {

    }
  }
}
