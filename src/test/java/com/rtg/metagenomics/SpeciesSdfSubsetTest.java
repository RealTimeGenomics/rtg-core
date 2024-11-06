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

package com.rtg.metagenomics;

import java.io.ByteArrayOutputStream;
import java.io.File;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.reader.Sdf2Fasta;
import com.rtg.taxonomy.TaxFilterCliTest;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

/**
 */
public class SpeciesSdfSubsetTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new SpeciesSdfSubset();
  }

  public void testHelp() {
    checkHelp("Extracts a subset of sequences from one SDF and outputs them to another SDF.",
        "output SDF",
        "input SDF",
        "file containing taxonomy");
  }

  public void testReduce() throws Exception {
    try (final TestDirectory dir = new TestDirectory("speciessdf")) {

      // create test metagenomic reference sdf
      final File fullSdf = TaxFilterCliTest.makeTaxonomySdf(dir);

      // call sdfsubset to create new sdf from original sdf and tax file
      final File reducedSdf = new File(dir, "sdf_reduced");
      final File reduceTaxonomy = new File(dir, "remove.tsv");
      FileHelper.resourceToFile("com/rtg/metagenomics/resources/taxonomy_remove.tsv", reduceTaxonomy);
      final String outStr = checkMainInitOk("-t", reduceTaxonomy.getAbsolutePath(), "-o", reducedSdf.getAbsolutePath(), "-i", fullSdf.getAbsolutePath());
      assertEquals("", outStr);

      // dump sequences to compare
      final Sdf2Fasta sdf2fasta = new Sdf2Fasta();

      final ByteArrayOutputStream out = new ByteArrayOutputStream();
      final MemoryPrintStream err = new MemoryPrintStream();
      final File reducedSequences = new File(dir, "sequences2.fasta");
      final int rc = sdf2fasta.mainInit(new String[] {"-i", reducedSdf.getAbsolutePath(), "-o", reducedSequences.getAbsolutePath(), "-Z"}, out, err.printStream());
      assertEquals("Error: " + err.toString(), "", err.toString());
      assertEquals(0, rc);
      //System.err.println(out.toString());

      // compare taxonomy tsv files in new sdf dir
      final File reducedTax = new File(reducedSdf, "taxonomy.tsv");
      final File reducedTaxLook = new File(reducedSdf, "taxonomy_lookup.tsv");

      mNano.check("taxonomy_remove.tsv", FileUtils.fileToString(reducedTax), false);
      mNano.check("taxonomy_lookup_reduced.tsv", FileUtils.fileToString(reducedTaxLook), false);
      mNano.check("sequences_reduced.fasta", FileUtils.fileToString(reducedSequences), false);
    }

  }

}
