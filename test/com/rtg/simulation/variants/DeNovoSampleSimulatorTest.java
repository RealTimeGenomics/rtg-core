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
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.reference.ReferenceGenome;
import com.rtg.reference.ReferenceGenome.ReferencePloidy;
import com.rtg.reference.Sex;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.PortableRandom;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.GenomePriorParamsBuilder;
import com.rtg.vcf.VcfUtils;

import junit.framework.TestCase;

/**
 */
public class DeNovoSampleSimulatorTest extends TestCase {

  @Override
  public void setUp() {
    Diagnostic.setLogStream();
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

  public void testSampleSimulator() throws IOException, InvalidParamsException {
    final File dir = FileUtils.createTempDir("genomemut2_", "test");
    try {
      final File sdf = new File(dir, "sdf");
      ReaderTestUtils.getDNADir(REF, sdf);
      final SequencesReader sr = SequencesReaderFactory.createMemorySequencesReader(sdf, true, LongRange.NONE);
      FileUtils.stringToFile(REFTXT, new File(sdf, ReferenceGenome.REFERENCE_FILE));

      // Generate variants
      final FixedStepPopulationVariantGenerator fixed = new FixedStepPopulationVariantGenerator(sr, 30, new Mutator("X"), new PortableRandom(10), 0.5);
      final List<PopulationVariantGenerator.PopulationVariant> variants = fixed.generatePopulation();
      final File popVcf = new File(dir, "popVcf.vcf.gz");
      PopulationVariantGenerator.writeAsVcf(popVcf, null, variants, sr);
      //String popVarStr = FileHelper.gzFileToString(popVcf);
      //System.out.println("-- Population Variants --");
      //System.out.println(popVarStr);

      // Generate a couple of samples w.r.t variants
      final SampleSimulator dadsim = new SampleSimulator(sr, new PortableRandom(42), ReferencePloidy.AUTO);
      final File dadVcf = new File(dir, "sample_dad.vcf.gz");
      dadsim.mutateIndividual(popVcf, dadVcf, "dad", Sex.MALE);
      final SampleSimulator momsim = new SampleSimulator(sr, new PortableRandom(43), ReferencePloidy.AUTO);
      final File momVcf = new File(dir, "sample_mom.vcf.gz");
      momsim.mutateIndividual(dadVcf, momVcf, "mom", Sex.FEMALE);

      final GenomePriorParams params = new GenomePriorParamsBuilder().create();

      // Now generate genotypes containing de novo variants
      final DeNovoSampleSimulator dad2sim = new DeNovoSampleSimulator(sr, params, new PortableRandom(63), ReferencePloidy.AUTO, 20, false);
      final File dad2Vcf = new File(dir, "sample_dad2.vcf.gz");
      dad2sim.mutateIndividual(momVcf, dad2Vcf, "dad", "dad2");

      final DeNovoSampleSimulator mom2sim = new DeNovoSampleSimulator(sr, params, new PortableRandom(64), ReferencePloidy.AUTO, 20, false);
      final File mom2Vcf = new File(dir, "sample_mom2.vcf.gz");
      mom2sim.mutateIndividual(dad2Vcf, mom2Vcf, "mom", "mom2");

      String sampleVcf = FileHelper.gzFileToString(mom2Vcf);
      //System.out.println("-- Final VCF --");
      //System.out.println(sampleVcf);
      sampleVcf = StringUtils.grepMinusV(sampleVcf, "^#");
      sampleVcf = StringUtils.grep(sampleVcf, VcfUtils.FORMAT_DENOVO);
      final String[] sampleVars = TestUtils.splitLines(sampleVcf);
      assertEquals(33, sampleVars.length); // Slightly less than 40 due to sex chromosomes, expect 35 +/- randomness.
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }
}
