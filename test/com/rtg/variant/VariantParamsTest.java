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
package com.rtg.variant;

import static com.rtg.util.StringUtils.LS;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import com.rtg.launcher.OutputParams;
import com.rtg.launcher.ReaderParams;
import com.rtg.launcher.SequenceParams;
import com.rtg.mode.SequenceMode;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reference.ReferenceGenome.ReferencePloidy;
import com.rtg.reference.Sex;
import com.rtg.sam.SharedSamConstants;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.LogRecord;
import com.rtg.util.io.LogStream;
import com.rtg.util.test.FileHelper;
import com.rtg.util.test.params.TestParams;

import junit.framework.TestCase;

/**
 */
public class VariantParamsTest extends TestCase {

  protected File mDir;

  @Override
  public void setUp() throws IOException {
    Diagnostic.setLogStream();
    mDir = FileHelper.createTempDirectory();
  }

  @Override
  public void tearDown() {
    Diagnostic.setLogStream();
    assertTrue(FileHelper.deleteAll(mDir));
    mDir = null;
  }

  private static final String TEST_OUTPUT = "varianttestoutput";


  private ReaderParams makeGenome() throws IOException {
    final File subjectsDir = FileUtils.createTempDir("test", "coverageparams", mDir);
    ReaderTestUtils.getReaderDNA(">t\nacgt", subjectsDir, null).close();
    return SequenceParams.builder().directory(subjectsDir).mode(SequenceMode.UNIDIRECTIONAL).create().readerParams();
  }

  public void testDefaultsWithNonIdentityPosterior() {
    final VariantParamsBuilder vpb = VariantParams.builder(); //.nonidentityPosterior(true);
    final VariantParams params = vpb.create();
    assertTrue(params.nonidentityPosterior());
    TestUtils.containsAll(params.toString(),
        "VariantParams q_default=20 mated_read_default=20 unmated_read_default=20" + LS,
        " mated_read_max=255 unmated_read_max=255 ignore_map_qualities=false" + LS,
        " hypercomplex_length=21 non_identity_posterior=true" + LS,
        " machine=null vcf_rp=false output_index=true" + LS,
        " call_level=INTERESTING indels=true no_complex_calls=false" + LS,
        " interesting_threshold=0.4343" + LS,
        " interesting_separation=4" + LS,
        " max_coverage_filter=2147483647:2147483647" + LS,
        " max_coverage_bypass=2147483647:2147483647" + LS,
        " ignore_quality_scores=false" + LS,
        " max_ambiguity=null" + LS,
        " sex=EITHER" + LS,
        " ploidy=AUTO" + LS,
        " chunk_size=1000 lookahead=2 max_read_length=1000" + LS,
        " threading_environment=PARALLEL treading_environment_seed=null" + LS,
        " exec_threads=1 io_threads=1" + LS,
        " hyper_complex_threshold=21" + LS,
        " ionTorrent=false prune_hypothesis=false trim_split=NONE" + LS,
        " somatic_rate=0.3 somatic_allele_balance=false somatic_ssp=none include_germline_variants=false include_gain_of_reference=false loh_prior=0.0 contamination_basis=1000" + LS,
        " no_disease_prior=0.95" + LS,
        " Relationships:null" + LS,
        " SamFilterParams minMapQ=-1 maxAlignmentCount=-1 maxMatedAlignmentScore=null maxUnmatedAlignmentScore=null",
        " AVR model: Not set" + LS,
        " min_avr_score=0.0" + LS,
        " genome_connectivity=null" + LS,
        " max_em_iterations=-1" + LS,
        " regions_bed_file=null" + LS
        );
  }

  public void testOutputParamsLink() throws Exception {
    final File outDir = new File(mDir, TEST_OUTPUT);
    assertTrue(outDir.mkdirs());
    final VariantParamsBuilder vpb = VariantParams.builder().outputParams(new OutputParams(outDir, false, false));
    VariantParams vp = vpb.create();
    assertEquals(outDir, vp.directory());
    assertEquals(new File(outDir, "test.txt"), vp.file("test.txt"));
    assertFalse(vp.blockCompressed());
    assertEquals(new File(outDir, "snps.vcf"), vp.vcfFile());
    vp.bedStream().close();
    assertTrue(vp.file("regions.bed").isFile());
    vp.vcfStream().close();
    assertTrue(vp.file("snps.vcf").isFile());
    vp = vpb.create();
    assertEquals(new File(outDir, "snps.vcf"), vp.vcfFile());
    vp.vcfStream().close();
    assertTrue(vp.file("snps.vcf").isFile());
  }

  public void testOk1() throws Exception {
    Diagnostic.setLogStream();
    final File map = File.createTempFile("testok1", "variantParams");
    try {
      FileUtils.stringToFile(SharedSamConstants.SAM9, map);
      assertTrue(map.isFile());

      VariantParams vp = null;
      try {
        final File outFile = new File(TEST_OUTPUT);
        assertTrue(outFile.mkdir());

        final String[] exp = {
            "VariantParams mapped reads=[" + map.getPath() + "]",
            " call_level=" + VariantOutputLevel.INTERESTING,
            " indels=" + Boolean.TRUE,
            " interesting_threshold=" + 5.0,
            " interesting_separation=" + 4,
            " heterozygous prior=0.0003000 homozygous prior=0.0003000",
            " SequenceParams mode=UNIDIRECTIONAL directory=",
            " OutputParams output directory=" + TEST_OUTPUT,
            " progress=" + Boolean.FALSE,
            " zip=" + Boolean.FALSE,
            " max_ambiguity=null",
            " sex=" + Sex.EITHER,
        };

        final List<File> mapped = new ArrayList<>();
        mapped.add(map);

        final VariantParamsBuilder vpb = VariantParams.builder();
        vpb.genomePriors(GenomePriorParams.builder().genomePriors("testhumanprior").create()).mapped(mapped).genome(makeGenome());
        vp = vpb.outputParams(new OutputParams(outFile, false, false)).interestingThreshold(5).create();

        final String ccs = vp.toString();
        //System.err.println(ccs);
        TestUtils.containsAll(ccs, exp);
        vp.close();
      } finally {
        if (vp != null) {
          FileHelper.deleteAll(vp.genome().directory());
        }
        FileHelper.deleteAll(new File(TEST_OUTPUT));
      }

    } finally {
      assertTrue(map.delete());
    }
  }

  //TODO do all calls
  public void testBuilder() throws InvalidParamsException, IOException {
    final VariantParamsBuilder vpb = VariantParams.builder();
    assertTrue(vpb.callLevel(VariantOutputLevel.ALL) == vpb);
    assertTrue(vpb.defaultQuality(42) == vpb);
    assertTrue(vpb.defaultMatedReadQuality(11) == vpb);
    assertTrue(vpb.defaultUnmatedReadQuality(12) == vpb);
    assertEquals(vpb.mapped(Collections.singleton(new File(""))), vpb);
    assertTrue(vpb.name("name") == vpb);
    assertTrue(vpb.interestingThreshold(5.0) == vpb);
    final OutputParams op = new OutputParams(new File(""), false, false);
    assertTrue(vpb.outputParams(op) == vpb);
    assertTrue(vpb.name("name") == vpb);
    assertTrue(vpb.genomePriors("human") == vpb);
    assertEquals(vpb.ioThreads(1), vpb);
    assertTrue(vpb.maxAmbiguity(1.0) == vpb);
    assertTrue(vpb.sex(Sex.FEMALE) == vpb);
    assertTrue(vpb.ploidy(ReferencePloidy.HAPLOID) == vpb);
    assertTrue(vpb.pruneHypotheses(true) == vpb);
    //    final SequenceParams sq = SequenceParams.builder().directory(new File("foo")).create();
    //    assertTrue(vpb.genome(null) == vpb);
    final VariantParams params = vpb.create();
    final String[] exp = {
        " q_default=42",
        " mated_read_default=11",
        " unmated_read_default=12",
        " call_level=" + VariantOutputLevel.ALL,
        " indels=" + Boolean.TRUE,
        " interesting_threshold=" + 5.0,
        " interesting_separation=" + 4,
        " heterozygous prior=0.0007100 homozygous prior=0.0005300",
        " progress=" + Boolean.FALSE,
        " zip=" + Boolean.FALSE,
        " max_ambiguity=" + Double.toString(1.0),
        " sex=" + Sex.FEMALE,
        " ploidy=" + ReferencePloidy.HAPLOID
    };
    final String actual = params.toString();
    //System.err.println(actual);
    TestUtils.containsAll(actual, exp);

    final VariantParams vpc = params.cloneBuilder().create();
    assertEquals(params.name(), vpc.name());
    assertEquals(op, vpc.outputParams());
    assertEquals(params.genome(), vpc.genome());
    assertEquals(params.mapped(), vpc.mapped());
    assertEquals(params.calibrations(), vpc.calibrations());
    assertEquals(params.qDefault(), vpc.qDefault());
    assertEquals(params.matedReadDefault(), vpc.matedReadDefault());
    assertEquals(params.unmatedReadDefault(), vpc.unmatedReadDefault());
    assertEquals(params.matedReadMax(), vpc.matedReadMax());
    assertEquals(params.unmatedReadMax(), vpc.unmatedReadMax());
    assertEquals(params.ignoreReadQualities(), vpc.ignoreReadQualities());
    assertEquals(params.filterParams(), vpc.filterParams());
    assertEquals(params.callLevel(), vpc.callLevel());
    assertEquals(params.outputNonSnps(), vpc.outputNonSnps());
    assertEquals(params.genomePriors(), vpc.genomePriors());
    assertEquals(params.interestingThreshold(), vpc.interestingThreshold());
    assertEquals(params.interestingSeparation(), vpc.interestingSeparation());
    assertEquals(params.hyperComplexLength(), vpc.hyperComplexLength());
    assertEquals(params.noComplexCalls(), vpc.noComplexCalls());
    assertEquals(params.maxCoverageFilter(), vpc.maxCoverageFilter());
    assertEquals(params.ignoreQualityScores(), vpc.ignoreQualityScores());
    assertEquals(params.ioThreads(), vpc.ioThreads());
    assertEquals(params.nonidentityPosterior(), vpc.nonidentityPosterior());
    assertEquals(params.maxAmbiguity(), vpc.maxAmbiguity());
    assertEquals(params.sex(), vpc.sex());
    assertEquals(params.machineErrorName(), vpc.machineErrorName());
    assertEquals(params.vcfRp(), vpc.vcfRp());
    assertEquals(params.outputIndex(), vpc.outputIndex());
    assertEquals(params.calibrator(), vpc.calibrator());
    assertEquals(params.chunkSize(), vpc.chunkSize());
    assertEquals(params.lookAhead(), vpc.lookAhead());
    assertEquals(params.ionTorrent(), vpc.ionTorrent());
    assertEquals(params.pruneHypotheses(), vpc.pruneHypotheses());
    assertEquals(params.ploidy(), vpc.ploidy());
    assertEquals(params.execThreads(), vpc.execThreads());
    assertEquals(params.maxEmIterations(), vpc.maxEmIterations());
    assertEquals(params.genomeConnectivity(), vpc.genomeConnectivity());
  }

  public void testStreams() throws IOException {
    final File fooFile = new File("foofile");
    try {
      final OutputParams op = new OutputParams(fooFile, false, false);
      final VariantParams params = VariantParams.builder().outputParams(op).create();

      final OutputStream os = params.vcfStream();
      assertNotNull(os);
      os.close();

      final OutputStream bs = params.bedStream();
      assertNotNull(bs);
      bs.close();
    } finally {
      FileUtils.deleteFiles(fooFile);
    }
  }

  /**
   * This tests the real human.properties values.
   * It will need to be updated each time they change.
   *
   * @throws InvalidParamsException
   * @throws IOException on IO error
   */
  public void testDefaults() throws InvalidParamsException, IOException {
    final VariantParamsBuilder vpb = VariantParams.builder();
    final VariantParams params = vpb.create();
    assertEquals(1, params.ioThreads());
    assertEquals(20, params.qDefault());
    assertEquals(20, params.matedReadDefault());
    assertEquals(20, params.unmatedReadDefault());
    assertEquals(255, params.matedReadMax());
    assertEquals(255, params.unmatedReadMax());
    assertEquals(Sex.EITHER, params.sex());
    assertEquals(ReferencePloidy.AUTO, params.ploidy());
    assertEquals("ModuleParams", params.name());
    assertTrue(params.nonidentityPosterior());
    assertFalse(params.ionTorrent());
    assertTrue(params.outputIndex());
    assertFalse(params.vcfRp());
    assertEquals(21, params.hyperComplexLength());
    assertFalse(params.ignoreReadQualities());
    assertEquals(0.0, params.somaticParams().lohPrior());
    assertEquals(-1, params.maxEmIterations());
    assertNull(params.genomeConnectivity());
  }

  /**
   * This tests that we can set values in the builder
   * and get the correct values back from the parameters object.
   *
   * @throws InvalidParamsException
   * @throws IOException on IO error
   */
  public void testNonDefaults() throws InvalidParamsException, IOException {
    final VariantParamsBuilder vpb = VariantParams.builder();
    vpb.ioThreads(4);
    vpb.defaultQuality(11);
    vpb.defaultMatedReadQuality(12);
    vpb.defaultUnmatedReadQuality(13);
    vpb.maxMatedReadQuality(42);
    vpb.maxUnmatedReadQuality(41);
    vpb.noComplexCalls(true);
    vpb.callLevel(VariantOutputLevel.ALL);
    vpb.interestingThreshold(3.5);
    vpb.interestingSeparation(4);
    vpb.chunkSize(5000);
    vpb.lookAhead(10);
    vpb.somaticParams(new SomaticParamsBuilder().lohPrior(0.1).create());
    vpb.maxEmIterations(42);
    vpb.genomeConnectivity(GenomeConnectivity.SPARSE);
    final VariantParams params = vpb.create();
    assertEquals("ModuleParams", params.name());
    assertEquals(4, params.ioThreads());
    assertEquals(11, params.qDefault());
    assertEquals(12, params.matedReadDefault());
    assertEquals(13, params.unmatedReadDefault());
    assertEquals(42, params.matedReadMax());
    assertEquals(41, params.unmatedReadMax());
    assertTrue(params.nonidentityPosterior());
    assertTrue(params.noComplexCalls());
    assertEquals(VariantOutputLevel.ALL, params.callLevel());
    assertEquals(3.5, params.interestingThreshold());
    assertEquals(4, params.interestingSeparation());
    assertEquals(5000, params.chunkSize());
    assertEquals(10, params.lookAhead());
    assertEquals(0.1, params.somaticParams().lohPrior());
    assertEquals(42, params.maxEmIterations());
    assertEquals(GenomeConnectivity.SPARSE, params.genomeConnectivity());
  }

  public void testBad1() throws InvalidParamsException, IOException {
    final LogStream logStream = new LogRecord();
    Diagnostic.setLogStream(logStream);
    final VariantParamsBuilder vpb = VariantParams.builder();
    try {
      vpb.genomePriors("badpriors");
      vpb.create();
      fail();
    } catch (final InvalidParamsException e) {
      assertEquals(ErrorType.INFO_ERROR, e.getErrorType());
      assertTrue(e.getMessage().contains("Invalid prior option \"badpriors\""));
    }
  }

  public void testBadChunkSize() {
    try {
      VariantParams.builder().chunkSize(0);
      fail();
    } catch (final IllegalArgumentException e) {
      // expected
    }
  }

  public void testBadLookAhead() {
    try {
      VariantParams.builder().lookAhead(1);
      fail();
    } catch (final IllegalArgumentException e) {
      // expected
    }
  }

  public void testMultisampleParams() {
    final VariantParams params = new VariantParamsBuilder().create();
    assertEquals(0.3, params.somaticParams().somaticRate(), 1e-8);
    assertEquals(0.95, params.noDiseasePrior(), 1e-8);
  }

  public void testOmnes() {
    new TestParams(VariantParams.class, VariantParamsBuilder.class).excludeParams("threadingLookAhead").check();
  }
}
