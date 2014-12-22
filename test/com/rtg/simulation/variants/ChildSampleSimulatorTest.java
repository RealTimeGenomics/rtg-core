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

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.List;

import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.reference.ReferenceGenome;
import com.rtg.reference.ReferenceGenome.DefaultFallback;
import com.rtg.reference.Sex;
import com.rtg.util.PortableRandom;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;
import com.rtg.util.test.NanoRegression;
import com.rtg.variant.util.MendeliannessChecker;

import junit.framework.TestCase;

/**
 */
public class ChildSampleSimulatorTest extends TestCase {

  NanoRegression mNano = null;
  @Override
  public void setUp() {
    Diagnostic.setLogStream();
    mNano = new NanoRegression(this.getClass());
  }

  @Override
  protected void tearDown() throws Exception {
    try {
      mNano.finish();
    } finally {
      mNano = null;
    }
  }

  private static final String REF =
      ">ref1" + StringUtils.LS
      + "cgtacattac" + "gagcgactag" + "ctagctagta" + "cgtacgtaca"
      + "atggcagcgt" + "attagcggca" + "aattgcgcat" + "tgcgtagcac"
      + "gcgcgattca" + "ttatgcgcgc" + "atcgatcgat" + "cgatcgatca" + StringUtils.LS
      + ">ref2" + StringUtils.LS
      + "atggcagcgt" + "attagcggca" + "aattgcgcat" + "tgcgtagcac"
      + "gcgcgattca" + "ttatgcgcgc" + "atcgatcgat" + "cgatcgatca"
      + "cgtacattac" + "gagcgactag" + "ctagctagta" + "cgtacgtaca" + StringUtils.LS
      + ">ref3" + StringUtils.LS
      + "gcgcgattca" + "ttatgcgcgc" + "atcgatcgat" + "cgatcgatca"
      + "cgtacattac" + "gagcgactag" + "ctagctagta" + "cgtacgtaca"
      + "atggcagcgt" + "attagcggca" + "aattgcgcat" + "tgcgtagcac" + StringUtils.LS
      + ">ref4" + StringUtils.LS
      + "cgtacattac" + "gagcgactag" + "ctagctagta" + "cgtacgtaca"
      + "gcgcgattca" + "ttatgcgcgc" + "atcgatcgat" + "cgatcgatca"
      + "atggcagcgt" + "attagcggca" + "aattgcgcat" + "tgcgtagcac" + StringUtils.LS;

  private static final String REFTXT =
      "# Simulated reference" + StringUtils.LS
      + "version 1" + StringUtils.LS
      + "either  def     diploid linear" + StringUtils.LS
      + "# ref1 ~= ChrX, ref2 ~= ChrY" + StringUtils.LS
      + "male    seq     ref1    haploid linear  ref2" + StringUtils.LS
      + "male    seq     ref2    haploid linear  ref1" + StringUtils.LS
      + "female  seq     ref1    diploid linear" + StringUtils.LS
      + "female  seq     ref2    none    linear" + StringUtils.LS
      + "either  seq     ref3    polyploid       circular" + StringUtils.LS;

  public void testSampleSimulator() throws IOException {
    try (final TestDirectory dir = new TestDirectory("childsim")) {
      final File sdf = new File(dir, "sdf");
      ReaderTestUtils.getDNADir(REF, sdf);
      final SequencesReader sr = SequencesReaderFactory.createMemorySequencesReader(sdf, true, LongRange.NONE);
      FileUtils.stringToFile(REFTXT, new File(sdf, ReferenceGenome.REFERENCE_FILE));

      // Generate variants
      FixedStepPopulationVariantGenerator fixed = new FixedStepPopulationVariantGenerator(sr, 10, new Mutator("X"), new PortableRandom(10), 0.5);
      List<PopulationVariantGenerator.PopulationVariant> variants = fixed.generatePopulation();
      final File popVcf = new File(dir, "popVcf.vcf.gz");
      sr.reset();
      PopulationVariantGenerator.writeAsVcf(popVcf, null, variants, sr);
      //String popVarStr = FileHelper.gzFileToString(popVcf);
      //System.out.println("-- Population Variants --");
      //System.out.println(popVarStr);

      // Generate sample w.r.t variants
      sr.reset();
      SampleSimulator dadsim = new SampleSimulator(sr, new PortableRandom(15), DefaultFallback.DIPLOID);
      File dadVcf = new File(dir, "sample_dad.vcf.gz");
      dadsim.mutateIndividual(popVcf, dadVcf, "dad", Sex.MALE);

      sr.reset();
      SampleSimulator momsim = new SampleSimulator(sr, new PortableRandom(65), DefaultFallback.DIPLOID);
      File momVcf = new File(dir, "sample_mom.vcf.gz");
      momsim.mutateIndividual(dadVcf, momVcf, "mom", Sex.FEMALE);

      // Generate children w.r.t variants
      sr.reset();
      ChildSampleSimulator sonsim = new ChildSampleSimulator(sr, new PortableRandom(76), DefaultFallback.DIPLOID, 0, false);
      File sonVcf = new File(dir, "sample_son.vcf.gz");
      sonsim.mutateIndividual(momVcf, sonVcf, "son", Sex.MALE, "dad", "mom");

      sr.reset();
      ChildSampleSimulator daughtersim = new ChildSampleSimulator(sr, new PortableRandom(13), DefaultFallback.DIPLOID, 0, false);
      File daughterVcf = new File(dir, "sample_daughter.vcf.gz");
      daughtersim.mutateIndividual(sonVcf, daughterVcf, "daughter", Sex.FEMALE, "dad", "mom");

      final ByteArrayOutputStream bos = new ByteArrayOutputStream();
      try {
        final MendeliannessChecker chk = new MendeliannessChecker();

        chk.mainInit(new String[] {"-t", sdf.getPath(), "-i", daughterVcf.getPath()}, bos, new PrintStream(bos));
      } finally {
        bos.close();
      }
      final String s = bos.toString().replaceAll("Checking: [^\n]*\n", "Checking: \n");
      TestUtils.containsAll(s,
          "Family: [dad + mom] -> [daughter, son]",
          "(0.00%) records did not conform to expected call ploidy",
          "(0.00%) records contained a violation of Mendelian constraints"
          );

      String sampleVcf = FileHelper.gzFileToString(daughterVcf);
      //System.out.println("-- Including sample foo --");
      //System.out.println(sampleVcf);
      sampleVcf = StringUtils.grepMinusV(sampleVcf, "^#");
      mNano.check("childsim", sampleVcf);
    }
  }
}
