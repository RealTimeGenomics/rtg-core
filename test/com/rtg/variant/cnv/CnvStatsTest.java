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
package com.rtg.variant.cnv;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.PrintStream;
import java.io.Reader;
import java.io.Writer;

import com.rtg.util.InvalidParamsException;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 * Tests corresponding class
 */
public class CnvStatsTest extends TestCase {
  private static final String TB = "\t";

  private File mDir = null;

  @Override
  public void setUp() throws IOException, InvalidParamsException {
    mDir = FileHelper.createTempDirectory();
    Diagnostic.setLogStream();
  }

  @Override
  public void tearDown() throws IOException  {
    assertTrue(!mDir.exists() || deleteAll(mDir));
    mDir = null;
  }

  static boolean deleteAll(final File file) {
    boolean ok = true;
    if (file != null) {
      if (file.isDirectory()) {
        final File[] files = file.listFiles();
        if (files != null) {
          for (final File x : files) {
            ok &= deleteAll(x);
          }
        }
      }
      ok &= file.delete();
    }
    return ok;
  }

  static final String HEADER = ""
    + "#Version" + StringUtils.LS
    + "#Seq" + TB + "position" + TB + "error" + TB + "bp-cn" + TB + "cn" + StringUtils.LS;

  private byte[] getCnvs(int[][] posLists, int[][] cnLists) {
    final ByteArrayOutputStream output = new ByteArrayOutputStream();
    final PrintStream cnvs = new PrintStream(output);
    final StringBuilder str = new StringBuilder();
    str.append(HEADER);
    for (int seqCount = 0; seqCount < posLists.length; ++seqCount) {
      final String seqName = "SEQ" + (seqCount + 1);
      final int[] positions = posLists[seqCount];
      final int[] cns = cnLists[seqCount];
      for (int j = 0; j < positions.length; ++j) {
        final int pos = positions[j];
        final int cn = cns[j];
        if (j + 1 < positions.length) {
          str.append(cnvLine(seqName, pos, positions[j + 1], cn));
        } else {
          str.append(cnvLine(seqName, pos, pos, cn));
        }
      }
    }
    cnvs.print(str.toString());
    cnvs.flush();
    cnvs.close();
    return output.toByteArray();
  }

  private Object cnvLine(String name, int pos, int pos2, int cn) {
    return name + TB + pos + TB + pos2 + TB + "cnv" + TB + cn + TB + "0.0" + TB + "0" + TB + StringUtils.LS;
  }

  /**
   * Test snp stats
   * @throws IOException if an IO Error occurs
   */
  public void testGetStats() throws IOException {
    //// tests one sequence
    // - test perfect
    int[] generatedPositions = {1, 100, 150, 200, 250, 300};
    int[] generatedCns = {2, 1, 3, 0, 1, 0};
    int[] detectedPositions = generatedPositions;
    int[] detectedCns = generatedCns;
    checkStats(new int[][]{generatedPositions}, new int[][]{generatedCns},
        new int[][]{detectedPositions}, new int[][]{detectedCns}, 6, 0, 0, false);

    // - test false positive
    generatedPositions = new int[] {1, 100, 150, 200, 250, 300};
    generatedCns = new int[] {2, 1, 3, 0, 1, 0};
    detectedPositions = new int[] {1, 100, 150, 200, 220, 250, 300};
    detectedCns = new int[] {2, 1, 3, 0, 2, 1, 0};
    checkStats(new int[][]{generatedPositions}, new int[][]{generatedCns},
        new int[][]{detectedPositions}, new int[][]{detectedCns}, 6, 1, 0, false);

    // - test false negative
    generatedPositions = new int[] {1, 100, 150, 200, 220, 250, 300};
    generatedCns = new int[] {2, 1, 3, 0, 2, 1, 0};
    detectedPositions = new int[] {1, 100, 150, 200, 250, 300};
    detectedCns = new int[] {2, 1, 3, 0, 1, 0};
    checkStats(new int[][]{generatedPositions}, new int[][]{generatedCns},
        new int[][]{detectedPositions}, new int[][]{detectedCns}, 6, 0, 1, false);

    // - test all false positive and false negative
    generatedPositions = new int[] {1, 100, 150, 200, 250, 300};
    generatedCns = new int[] {2, 1, 3, 0, 1, 0};
    detectedPositions = new int[] {1, 70, 120, 170, 220, 270, 300};
    detectedCns = new int[] {2, 2, 1, 3, 0, 1, 0};
    checkStats(new int[][]{generatedPositions}, new int[][]{generatedCns},
        new int[][]{detectedPositions}, new int[][]{detectedCns}, 2, 5, 4, false);
  }

  public void testGetStatsMulti() throws IOException {
    //// tests multiple sequences
    // - test perfect
    final int[][] generatedSequences = {{1, 100, 150, 200, 250, 300},
        {1, 600, 2001}};

    final int[][] generatedSequenceCNs = {{2, 1, 3, 0, 1, 4},
        {1, 7, 0}};
    checkStats(generatedSequences, generatedSequenceCNs
        , generatedSequences, generatedSequenceCNs, 9, 0, 0, false);
  }

  public void testGetStatsMultiCollapse() throws IOException {
    //// tests multiple sequences
    // - test perfect
    final int[][] generatedSequences = {{1, 100, 150, 200, 250, 300},
        {1, 600, 2001}};

    final int[][] generatedSequenceCNs = {{2, 1, 3, 3, 1, 4},
        {1, 7, 0}};
    checkStats(generatedSequences, generatedSequenceCNs, generatedSequences, generatedSequenceCNs, 8, 1, 0, true);
  }

  public void checkStats(int[][] generatedPositions, int[][] generatedCns,
      int[][] detectedPositions, int[][] detectedCns,
      final int correctPos, final int falsePositive, final int falseNegative, boolean collapse) throws IOException {
    final byte[] generatedInput = getCnvs(generatedPositions, generatedCns);
    final byte[] detectedInput = getCnvs(detectedPositions, detectedCns);

    final Reader generated = new InputStreamReader(new ByteArrayInputStream(generatedInput));
    final Reader detected = new InputStreamReader(new ByteArrayInputStream(detectedInput));
    final File outFile = new File(mDir, "out");
    final Writer outStream = new OutputStreamWriter(new FileOutputStream(new File(mDir, "out")));
    try {
      final CnvStats stats = new CnvStats();
      stats.setCollapse(collapse);
      stats.setThreshold(10);
      stats.getStats(generated, detected, outStream);
      outStream.flush();
      final String outString = FileUtils.fileToString(outFile);
      checkResults(outString, correctPos, falsePositive, falseNegative);
    } finally {
      outStream.flush();
      outStream.close();
    }
  }

  private void checkResults(String outString, int correctPos, int falsePositive, int falseNegative) {
    final String[] lines = outString.split(StringUtils.LS);
    //System.err.println("num lines " + lines.length + "; num sequences " + numSequences);
    final String[] parts = lines[2].split(TB);
    assertEquals("correct-pos", correctPos, Integer.parseInt(parts[0]));
    assertEquals("false-positive", falsePositive, Integer.parseInt(parts[1]));
    assertEquals("false-negative", falseNegative, Integer.parseInt(parts[2]));
  }

}
