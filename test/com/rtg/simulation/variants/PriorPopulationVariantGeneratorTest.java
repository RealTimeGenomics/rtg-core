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
package com.rtg.simulation.variants;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.PortableRandom;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;
import com.rtg.util.test.NanoRegression;
import com.rtg.variant.GenomePriorParams;

import junit.framework.TestCase;

/**
 */
public class PriorPopulationVariantGeneratorTest extends TestCase {

  private NanoRegression mNano = null;

  @Override
  public void setUp() {
    Diagnostic.setLogStream();
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

  private static final String REF = ">ref" + StringUtils.LS
          + "cgtacattac" + "gagcgactag" + "ctagctagta" + "cgtacgtaca"
          + "atggcagcgt" + "attagcggca" + "aattgcgcat" + "tgcgtagcac"
          + "gcgcgattca" + "ttatgcgcgc" + "atcgatcgat" + "cgatcgatca"
          + "cgtacattac" + "gagcgactag" + "ctagctagta" + "cgtacgtaca"
          + "atggcagcgt" + "attagcggca" + "aattgcgcat" + "tgcgtagcac"
          + "atggcagcgt" + "attagcggca" + "aattgcgcat" + "tgcgtagcac"
          + "gcgcgattca" + "ttatgcgcgc" + "atcgatcgat" + "cgatcgatca"
          + "cgtacattac" + "gagcgactag" + "ctagctagta" + "cgtacgtaca"
          + "atggcagcgt" + "attagcggca" + "aattgcgcat" + "tgcgtagcac"
          + "gcgcgattca" + "ttatgcgcgc" + "atcgatcgat" + "cgatcgatca";

  public void testVariantGenerator() throws IOException, InvalidParamsException {
    final File dir = FileUtils.createTempDir("genomemut2_", "test");
    try {
      final SequencesReader sr = ReaderTestUtils.getReaderDnaMemory(REF);

      // Make some priors that will let things happen more often
      final GenomePriorParams priors = GenomePriorParams.builder()
      .genomeSnpRateHetero(0.05).genomeSnpRateHomo(0.05)           // The priors probably should have just one rather than splitting
      .genomeMnpBaseRateHetero(0.05).genomeMnpBaseRateHomo(0.05)   // The priors probably should have just one rather than splitting
      .genomeIndelEventRate(0.05)
      .create();

      // Generate variants
      final PriorPopulationVariantGenerator gen = new PriorPopulationVariantGenerator(sr, new PopulationMutatorPriors(priors), new PortableRandom(10), 1);
      final List<PopulationVariantGenerator.PopulationVariant> variants = gen.generatePopulation();
      final File popVcf = new File(dir, "popVcf.vcf.gz");
      PopulationVariantGenerator.writeAsVcf(popVcf, null, variants, sr);
      final String popVarStr = FileHelper.gzFileToString(popVcf);
      //System.out.println("-- Population Variants --");
      //System.out.println(popVarStr);
      final String vars = StringUtils.grepMinusV(popVarStr, "^#");
      mNano.check("popvars", vars);
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }
}
