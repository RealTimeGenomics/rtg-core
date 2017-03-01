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
import java.util.Arrays;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.launcher.MainResult;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

/**
 */
public class TaxFilterCliTest extends AbstractCliTest {

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

  public File makeTaxonomy(File dir) throws IOException {
    final File taxonomy = new File(dir, "tax.tsv");
    FileHelper.resourceToFile("com/rtg/taxonomy/resources/taxonomy.tsv", taxonomy);
    return taxonomy;
  }

  public static File makeTaxonomySdf(File dir) throws IOException {
    final File fullSdf = new File(dir, "sdf_full");

    // Format the sequence data
    ReaderTestUtils.getDNADir(FileHelper.resourceToString("com/rtg/metagenomics/resources/sequences.fasta"), fullSdf);

    // cp taxonomy files to sdf
    FileHelper.resourceToFile("com/rtg/metagenomics/resources/taxonomy.tsv", new File(fullSdf, "taxonomy.tsv"));
    FileHelper.resourceToFile("com/rtg/metagenomics/resources/taxonomy_lookup.tsv", new File(fullSdf, "taxonomy_lookup.tsv"));
    return fullSdf;
  }

  public MainResult checkResults(TestDirectory dir, String id, boolean sdf, String... extraargs) throws IOException {
    final File output = new File(dir, id + (sdf ? "_sdf" : "_tsv") + "_output");
    final File taxInput = sdf ? makeTaxonomySdf(dir) : makeTaxonomy(dir);
    String[] args = {
      "-i", taxInput.getAbsolutePath(), "-o", output.getPath()
    };
    args = Arrays.copyOf(args, args.length + extraargs.length);
    System.arraycopy(extraargs, 0, args, args.length - extraargs.length, extraargs.length);

    final MainResult res = checkMainInit(args);

    mNano.check("taxonomy_" + id + ".tsv", FileUtils.fileToString(sdf ? new File(output, TaxonomyUtils.TAXONOMY_FILE) : output), false);
    if (sdf) {
      mNano.check("taxonomy_" + id + "_lookup.tsv", FileUtils.fileToString(new File(output, TaxonomyUtils.TAXONOMY_TO_SEQUENCE_FILE)), false);
    }

    return res;
  }

  public void testSubset() throws Exception {
    try (TestDirectory dir = new TestDirectory()) {
      final File subsetids = new File(dir, "subset.txt");
      FileHelper.resourceToFile("com/rtg/taxonomy/resources/ids_subset.txt", subsetids);
      checkResults(dir, "subset", false, "-s", subsetids.getPath());
      checkResults(dir, "subset", true, "-s", subsetids.getPath());
    }
  }

  public void testSubTree() throws Exception {
    try (TestDirectory dir = new TestDirectory()) {
      final File subtreeids = new File(dir, "subtree.txt");
      FileHelper.resourceToFile("com/rtg/taxonomy/resources/ids_subtree.txt", subtreeids);
      checkResults(dir, "subtree", false, "-S", subtreeids.getPath());
      checkResults(dir, "subtree", true, "-S", subtreeids.getPath());
    }
  }

  public void testRemove() throws Exception {
    try (TestDirectory dir = new TestDirectory()) {
      final File removeids = new File(dir, "remove.txt");
      FileHelper.resourceToFile("com/rtg/taxonomy/resources/ids_remove.txt", removeids);

      checkResults(dir, "remove", false, "-r", removeids.getPath());
      checkResults(dir, "remove", true, "-r", removeids.getPath());
    }
  }

  public void testRemoveSequences() throws Exception {
    try (TestDirectory dir = new TestDirectory()) {
      final File removeids = new File(dir, "remove.txt");
      FileHelper.resourceToFile("com/rtg/taxonomy/resources/ids_remove_seqs.txt", removeids);

      checkResults(dir, "remove_sequences", true, "-R", removeids.getPath());

      FileHelper.resourceToFile("com/rtg/taxonomy/resources/ids_remove_all_seqs.txt", removeids);
      checkResults(dir, "remove_all_sequences", true, "-R", removeids.getPath());
    }
  }

  public void testRenameNoRank() throws Exception {
    try (TestDirectory dir = new TestDirectory()) {
      final File rename = new File(dir, "remove.txt");
      FileHelper.resourceToFile("com/rtg/taxonomy/resources/no_rank_rename.txt", rename);

      final MainResult res = checkResults(dir, "norank_rename", false, "--rename-norank", rename.getPath());

      TestUtils.containsAll(res.err(),
        "Node not found in taxonomy: 999999",
        "Node 1118 rank is not \"no rank\": order");
    }
  }

}
