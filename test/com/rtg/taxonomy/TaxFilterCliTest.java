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

package com.rtg.taxonomy;

import java.io.File;
import java.io.IOException;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;
import com.rtg.util.test.NanoRegression;

/**
 */
public class TaxFilterCliTest extends AbstractCliTest {
  private NanoRegression mNano;

  @Override
  public void setUp() throws IOException {
    super.setUp();
    mNano = new NanoRegression(TaxFilterCliTest.class);
  }

  @Override
  public void tearDown() throws IOException {
    super.tearDown();
    try {
      mNano.finish();
    } finally {
      mNano = null;
    }
  }

  @Override
  public AbstractCli getCli() {
    return new TaxFilterCli();
  }

  public void testHelp() {
    checkHelp("Reference taxonomy filtering.",
              "output=FILE",
              "input=FILE",
              "subset=FILE",
              "remove=FILE",
              "rename-norank=FILE",
              "taxonomy input");
  }

  public void testSubset() throws Exception {
    try (TestDirectory dir = new TestDirectory()) {
      final File taxonomy = new File(dir, "tax.tsv");
      FileHelper.resourceToFile("com/rtg/taxonomy/resources/taxonomy.tsv", taxonomy);
      final File subsetids = new File(dir, "subset.txt");
      FileHelper.resourceToFile("com/rtg/taxonomy/resources/ids_subset.txt", subsetids);
      final File output = new File(dir, "output.txt");

      final String out = checkMainInitOk("-i", taxonomy.getAbsolutePath(), "-o", output.getPath(), "-s", subsetids.getPath());
      assertEquals("", out);
      final String actual = FileUtils.fileToString(output);
      mNano.check("taxonomy_subset.tsv", actual, false);
    }
  }

  public void testRemove() throws Exception {
    try (TestDirectory dir = new TestDirectory()) {
      final File taxonomy = new File(dir, "tax.tsv");
      FileHelper.resourceToFile("com/rtg/taxonomy/resources/taxonomy.tsv", taxonomy);
      final File removeids = new File(dir, "remove.txt");
      FileHelper.resourceToFile("com/rtg/taxonomy/resources/ids_remove.txt", removeids);
      final File output = new File(dir, "output.txt");

      final String out = checkMainInitOk("-i", taxonomy.getAbsolutePath(), "-o", output.getPath(), "-r", removeids.getPath());
      assertEquals("", out);
      final String actual = FileUtils.fileToString(output);
      mNano.check("taxonomy_remove.tsv", actual, false);
    }
  }

  public void testRenameNoRank() throws Exception {
    try (TestDirectory dir = new TestDirectory()) {
      final File taxonomy = new File(dir, "tax.tsv");
      FileHelper.resourceToFile("com/rtg/taxonomy/resources/taxonomy.tsv", taxonomy);
      final File rename = new File(dir, "remove.txt");
      FileHelper.resourceToFile("com/rtg/taxonomy/resources/no_rank_rename.txt", rename);
      final File output = new File(dir, "output.txt");

      final String out = checkMainInitWarn("-i", taxonomy.getAbsolutePath(), "-o", output.getPath(), "--rename-norank", rename.getPath());
      //assertEquals("", out);
      assertTrue(out.contains("Node not found in taxonomy: 999999"));
      assertTrue(out.contains("Node 1118 rank is not \"no rank\": order"));
      final String actual = FileUtils.fileToString(output);
      mNano.check("taxonomy_norank_rename.tsv", actual, false);
    }
  }

}
