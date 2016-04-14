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

package com.rtg.blacklist;

import java.io.File;
import java.io.IOException;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.launcher.MainResult;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.NanoRegression;

/**
 *
 */
public class HashDistCliTest extends AbstractCliTest {

  private NanoRegression mNano;

  @Override
  public void setUp() throws IOException {
    super.setUp();
    mNano = new NanoRegression(HashDistCliTest.class);
  }

  @Override
  public void tearDown() throws IOException {
    super.tearDown();
    mNano.finish();
    mNano = null;
  }

  @Override
  protected AbstractCli getCli() {
    return new HashDistCli();
  }

  private static final String REF = ">1\nacagctaatgcat"; //a 5 c 3 g 2 t 3

  public void test() throws IOException {
    try (TestDirectory dir = new TestDirectory()) {
      final File sdfDir = ReaderTestUtils.getDNADir(REF, new File(dir, "sdf"));
      final File outDir = new File(dir, "out");
      final MainResult mr = MainResult.run(getCli(), sdfDir.getPath(), "--word", "1", "--output", outDir.getPath(), "--blacklist-threshold", "1");
      final String histogram = FileUtils.fileToString(new File(outDir, "histogram.txt"));
      final String blacklist = FileUtils.fileToString(new File(outDir, "blacklist"));
      mNano.check("histogram.txt", histogram);
      mNano.check("blacklist.txt", blacklist);
    }
  }
}