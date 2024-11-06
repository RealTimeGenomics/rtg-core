/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
