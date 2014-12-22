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

package com.rtg.metagenomics;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.launcher.GlobalFlags;
import com.rtg.reader.FormatCli;
import com.rtg.reader.Sdf2Fasta;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.test.FileHelper;
import com.rtg.util.test.NanoRegression;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 */
public class SpeciesSdfSubsetTest extends AbstractCliTest {
  private NanoRegression mNano;

  @Override
  public void setUp() throws IOException {
    super.setUp();
    GlobalFlags.resetAccessedStatus();
    mNano = new NanoRegression(SpeciesSdfSubsetTest.class);
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

  public static Test suite() {
    return new TestSuite(SpeciesSdfSubsetTest.class);
  }

  public static void main(final String[] args) {
    junit.textui.TestRunner.runAndWait(suite());
  }

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
    final File dir = FileUtils.createTempDir("speciessdf", "test");
    try {

      // create full sdf from sequences.fasta
      final File fullSdf = new File(dir, "sdf_full");
      final File sequences = new File(dir, "sequences.fasta");
      FileHelper.resourceToFile("com/rtg/metagenomics/resources/sequences.fasta", sequences);
      ByteArrayOutputStream out = new ByteArrayOutputStream();
      MemoryPrintStream err = new MemoryPrintStream();
      final FormatCli format = new FormatCli();
      int rc = format.mainInit(new String[] {"-o", fullSdf.getAbsolutePath(), sequences.getAbsolutePath()}, out, err.printStream());
      assertEquals("Error: " + err.toString(), "", err.toString());
      assertEquals(0, rc);
      //System.err.println(out.toString());

      // cp taxonomy files to sdf

      final File taxonomy = new File(fullSdf, "taxonomy.tsv");
      FileHelper.resourceToFile("com/rtg/metagenomics/resources/taxonomy.tsv", taxonomy);
      final File taxonomyLookup = new File(fullSdf, "taxonomy_lookup.tsv");
      FileHelper.resourceToFile("com/rtg/metagenomics/resources/taxonomy_lookup.tsv", taxonomyLookup);

      // call sdfsubset to create new sdf from original sdf and tax file
      final File reducedSdf = new File(dir, "sdf_reduced");
      final File reduceTaxonomy = new File(dir, "remove.tsv");
      FileHelper.resourceToFile("com/rtg/metagenomics/resources/taxonomy_remove.tsv", reduceTaxonomy);
      final String outStr = checkMainInitOk("-t", reduceTaxonomy.getAbsolutePath(), "-o", reducedSdf.getAbsolutePath(), "-i", fullSdf.getAbsolutePath());
      assertEquals("", outStr);

      // dump sequences to compare
      final Sdf2Fasta sdf2fasta = new Sdf2Fasta();

      out = new ByteArrayOutputStream();
      err = new MemoryPrintStream();
      final File reducedSequences = new File(dir, "sequences2.fasta");
      rc = sdf2fasta.mainInit(new String[] {"-i", reducedSdf.getAbsolutePath(), "-o", reducedSequences.getAbsolutePath(), "-Z"}, out, err.printStream());
      assertEquals("Error: " + err.toString(), "", err.toString());
      assertEquals(0, rc);
      //System.err.println(out.toString());

      // compare taxonomy tsv files in new sdf dir
      final File reducedTax = new File(reducedSdf, "taxonomy.tsv");
      final File reducedTaxLook = new File(reducedSdf, "taxonomy_lookup.tsv");

      mNano.check("taxonomy_remove.tsv", FileUtils.fileToString(reducedTax), false);
      mNano.check("taxonomy_lookup_reduced.tsv", FileUtils.fileToString(reducedTaxLook), false);
      mNano.check("sequences_reduced.fasta", FileUtils.fileToString(reducedSequences), false);
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }

  }

}
