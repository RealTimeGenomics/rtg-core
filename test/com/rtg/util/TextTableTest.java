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

import java.io.IOException;

import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.SlimException;
import com.rtg.util.test.NanoRegression;

import junit.framework.TestCase;

/**
 */
public class TextTableTest extends TestCase {

  NanoRegression mNano;

  @Override
  public void setUp() throws IOException {
    mNano = new NanoRegression(this.getClass());
  }

  @Override
  public void tearDown() throws IOException {
    try {
      mNano.finish();
    } finally {
      mNano = null;
    }
  }


  public void testDefault() throws IOException {
    final TextTable table = new TextTable();
    assertEquals("", table.toString());
    assertEquals(0, table.numRows());

    table.addRow("depth", "breadth", "covered", "size", "non-N-depth", "non-N-breadth", "non-N-covered", "non-N-size", "name");
    table.addSeparator();
    table.addRow("0.0097", "0.0087", "2140198", "247249719", "0.0107", "0.0095", "2140198", "224999719", "chr1");
    table.addRow("0.0112", "0.0093", "1259191", "135374737", "0.0116", "0.0096", "1259191", "131624728", "chr10");
    table.addRow("0.0100", "0.0093", "1256531", "134452384", "0.0102", "0.0096", "1256531", "131130753", "chr11");
    table.addRow("0.0093", "0.0088", "1162746", "132349534", "0.0095", "0.0089", "1162746", "130303032", "chr12");
    assertEquals(6, table.numRows());

    mNano.check("tt-default.txt", table.toString());

    table.setAlignment(TextTable.Align.LEFT, TextTable.Align.RIGHT, TextTable.Align.CENTER, TextTable.Align.LEFT, TextTable.Align.CENTER, TextTable.Align.RIGHT, TextTable.Align.LEFT);
    mNano.check("tt-mixed-alignments.txt", table.toString());
  }

  public void testIndent() throws IOException {
    final TextTable table = new TextTable(3, 2, TextTable.Align.RIGHT);
    table.addRow("1:", "4");
    table.addRow("2:", "42");
    table.addRow("10:", "24");
    table.addRow("longer:", "333");

    mNano.check("tt-indent.txt", table.toString());

    table.setAlignment(TextTable.Align.LEFT);
    mNano.check("tt-indent-1left.txt", table.toString());

    table.setAlignment(TextTable.Align.CENTER);
    mNano.check("tt-indent-1center.txt", table.toString());
  }

  public void testSingleRowAndColumn() {
    final TextTable table = new TextTable();
    table.addRow("1");
    assertEquals("1" + StringUtils.LS, table.toString());
  }

  public void testErrors() {
    Diagnostic.setLogStream();
    final TextTable table = new TextTable();
    table.addRow("1", "2");
    try {
      table.addRow("1", "2", "3");
      fail();
    } catch (SlimException e) {
      assertEquals("Mismatching number of columns in table formatter", e.getMessage());
    }
    try {
      table.addRow("1", null);
      fail();
    } catch (SlimException e) {
      assertEquals("Null provided as value in table formatters", e.getMessage());
    }
  }
}
