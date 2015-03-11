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
package com.rtg.variant.bayes.multisample;

import static com.rtg.sam.SharedSamConstants.OUT_SAM;
import static com.rtg.sam.SharedSamConstants.REF_SEQS;
import static com.rtg.sam.SharedSamConstants.REF_SEQS11;
import static com.rtg.sam.SharedSamConstants.SAM11;
import static com.rtg.sam.SharedSamConstants.SAM_FAMILY;
import static com.rtg.util.StringUtils.FS;
import static com.rtg.util.StringUtils.LS;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;

import com.rtg.launcher.AbstractNanoTest;
import com.rtg.launcher.MainResult;
import com.rtg.launcher.OutputParams;
import com.rtg.launcher.SequenceParams;
import com.rtg.mode.SequenceMode;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reference.Sex;
import com.rtg.relation.GenomeRelationships;
import com.rtg.relation.Relationship.RelationshipType;
import com.rtg.sam.SamUtils;
import com.rtg.sam.SharedSamConstants;
import com.rtg.tabix.TabixIndexer;
import com.rtg.tabix.UnindexableDataException;
import com.rtg.usage.UsageMetric;
import com.rtg.util.HtmlReportHelper;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.TestUtils;
import com.rtg.util.Utils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.BgzipFileHelper;
import com.rtg.util.test.FileHelper;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.StaticThreshold;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantParamsBuilder;
import com.rtg.variant.bayes.multisample.population.PopulationCallerConfiguration;
import com.rtg.variant.bayes.multisample.population.PopulationNanoTest;
import com.rtg.variant.bayes.multisample.singleton.SingletonCli;
import com.rtg.vcf.VariantStatistics;

/**
 */
public class MultisampleTaskTest extends AbstractNanoTest {

  private static final String RESOURCE_DIR = "com/rtg/variant/bayes/multisample/resources/";

  protected VariantParamsBuilder getBuilder() {
    return new VariantParamsBuilder();
  }

  protected String[] getArguments() {
    return new String[0];
  }

  public void testCancer() throws IOException, InvalidParamsException, UnindexableDataException {
    final GenomeRelationships genomeRelationships = new GenomeRelationships();
    genomeRelationships.addGenome("cancer");
    genomeRelationships.addGenome("TEST");
    genomeRelationships.addRelationship(RelationshipType.ORIGINAL_DERIVED, "TEST", "cancer").setProperty("contamination", "0.3");
    check(SharedSamConstants.SAM_CANCER, genomeRelationships, 6 * 8L);
  }

  public void testFamily() throws IOException, InvalidParamsException, UnindexableDataException {
    final GenomeRelationships genomeRelationships = new GenomeRelationships();
    genomeRelationships.addGenome("father", Sex.MALE);
    genomeRelationships.addGenome("mother", Sex.FEMALE);
    genomeRelationships.addGenome("child", Sex.MALE);
    genomeRelationships.addParentChild("father", "child");
    genomeRelationships.addParentChild("mother", "child");
    check(SAM_FAMILY, genomeRelationships, 6 * 8L);
  }

  private void check(final String sam, final GenomeRelationships genomeRelationships, long usageExp) throws IOException, UnindexableDataException, InvalidParamsException {
    final VariantParamsBuilder builder = getBuilder();
    check(sam, genomeRelationships, builder, usageExp);
  }

  private void check(final String sam, final GenomeRelationships genomeRelationships, final VariantParamsBuilder builder, long usageExp) throws IOException, UnindexableDataException, InvalidParamsException {
    Diagnostic.setLogStream();
    try (final TestDirectory input = new TestDirectory("variance_in");
         final TestDirectory outFile = new TestDirectory("variance_out")) {
      final File samFile = new File(input, "alignments.sam.gz");
      FileHelper.stringToGzFile(sam, samFile);
      new TabixIndexer(samFile, new File(input, "alignments.sam.gz.tbi")).saveSamIndex();
      final File templ = ReaderTestUtils.getDNADir(REF_SEQS, FileUtils.createTempDir("testfasta", "dna"), true, true, true);
      final SequenceParams g = SequenceParams.builder().directory(templ).mode(SequenceMode.UNIDIRECTIONAL).create();
      try (ByteArrayOutputStream bos = new ByteArrayOutputStream()) {
        final ArrayList<File> mapped = new ArrayList<>();
        mapped.add(samFile);
        builder.mapped(mapped);
        builder.genome(g.readerParams());
        builder.outputParams(new OutputParams(outFile, false, false));
        builder.genomePriors(GenomePriorParams.builder().contraryProbability(1).create());
        builder.genomeRelationships(genomeRelationships);
        builder.maxCoverageFilter(new StaticThreshold(100));
        builder.uberHeader(SamUtils.getUberHeader(null, mapped, false, genomeRelationships == null ? null : genomeRelationships.genomes()));
        final VariantParams p = builder.create();
        final UsageMetric usageMetric = new UsageMetric();
        final MultisampleTask<VariantStatistics> task = new MultisampleTask<>(p, new PopulationCallerConfiguration.Configurator(), bos, new VariantStatistics(null), usageMetric);
        final ByteArrayOutputStream berr = new ByteArrayOutputStream();
        try {
          final PrintStream err = new PrintStream(berr);
          try {
            //Diagnostic.setLogStream(System.err);
            task.exec();
          } finally {
            err.flush();
            err.close();
          }
          final String errStr = berr.toString();
          assertEquals(errStr, "", errStr);
          assertEquals(usageExp, task.usage());
        } finally {
          FileHelper.deleteAll(templ);
        }
      }
      // final String s = bos.toString();
      // System.out.println(s);
    }
  }

  private static final String EXP_STAT = ""
      + "Failed Filters               : 0" + LS
      + "Passed Filters               : 2" + LS
      + "SNPs                         : 2" + LS
      + "MNPs                         : 0" + LS
      + "Insertions                   : 0" + LS
      + "Deletions                    : 0" + LS
      + "Indels                       : 0" + LS
      + "Same as reference            : 0" + LS
      + "SNP Transitions/Transversions: 0.00 (0/4)" + LS
      + "Total Het/Hom ratio          : - (2/0)" + LS
      + "SNP Het/Hom ratio            : - (2/0)" + LS
      + "MNP Het/Hom ratio            : - (0/0)" + LS
      + "Insertion Het/Hom ratio      : - (0/0)" + LS
      + "Deletion Het/Hom ratio       : - (0/0)" + LS
      + "Indel Het/Hom ratio          : - (0/0)" + LS
      + "Insertion/Deletion ratio     : - (0/0)" + LS
      + "Indel/SNP+MNP ratio          : 0.00 (0/2)" + LS;

  public void testStats() throws Exception {
    //REF_SEQS11, SAM11, args0, 11, "", 0
    //Diagnostic.setLogStream(System.err);
    try (final TestDirectory outroot = new TestDirectory("variance_out");
         final TestDirectory input = new TestDirectory("variance_in")) {
      // FileUtils.saveFile(new File(input, OUT_SAM), sam);
      final File file = new File(input, OUT_SAM + ".gz");
      final File output = new File(outroot, "out");
      BgzipFileHelper.streamToBgzipFile(new ByteArrayInputStream(SAM11.getBytes()), file);
      new TabixIndexer(file, new File(input, OUT_SAM + ".gz.tbi")).saveSamIndex();
      final File templ = ReaderTestUtils.getDNADir(REF_SEQS11, new File(input, "dna"));
      final String[] args = {
        "-t", templ.getPath(),
        "-o", output.getPath(),
        "-Z",
        "-m", "default", "--keep-duplicates",
        input.getPath() + FS + OUT_SAM + ".gz",
        "--" + AbstractMultisampleCli.NO_CALIBRATION
      };
      final MainResult res = MainResult.run(new SingletonCli(), Utils.append(args, getArguments()));
      assertEquals("", res.err());
      assertEquals(EXP_STAT, FileUtils.fileToString(new File(output, "summary.txt")));

      // Test some aspects of html report
      final HtmlReportHelper helper = new HtmlReportHelper(output, "index");
      final String report = FileUtils.fileToString(helper.getReportFile());
      assertFalse(report.contains("Sample:"));
      assertTrue(report.contains("0.00 (0/4)"));

      final String log = FileUtils.fileToString(new File(output, "snp.log"));
      //System.err.println(log);
      TestUtils.containsAll(log, "Timer INCR");
    }
  }

  public void testNonExistentRef() throws Exception {
    //REF_SEQS11, SAM11, args0, 11, "", 0
    //Diagnostic.setLogStream(System.err);
    try (final TestDirectory outroot = new TestDirectory("variance_out");
         final TestDirectory input = new TestDirectory("variance_in")) {
      // FileUtils.saveFile(new File(input, OUT_SAM), sam);
      final File file = new File(input, OUT_SAM + ".gz");
      final File output = new File(outroot, "out");
      BgzipFileHelper.streamToBgzipFile(new ByteArrayInputStream(SAM11.getBytes()), file);
      new TabixIndexer(file, new File(input, OUT_SAM + ".gz.tbi")).saveSamIndex();
      final File templ = ReaderTestUtils.getDNADir(REF_SEQS11, new File(input, "dna"));
      final String[] args = {
        "-t", templ.getPath(),
        "-o", output.getPath(),
        "-Z",
        "-m", "default", "--keep-duplicates",
        input.getPath() + FS + OUT_SAM + ".gz",
        "--region", "foooo",
        "--" + AbstractMultisampleCli.NO_CALIBRATION
      };
      final MainResult res = MainResult.run(new SingletonCli(), Utils.append(args, getArguments()));
      assertEquals(1, res.rc());
      assertEquals("Error: Sequence \"foooo\" referenced in region was not found in the SAM sequence dictionary." + LS, res.err());
    }
  }

  // Leverage the population caller nano regression code to test indel triggering
  public void testIndelTriggering() throws Exception {
    final String testName = "indel_triggering_fraction";
    final String sam = FileHelper.resourceToString(RESOURCE_DIR + testName + ".sam");
    PopulationNanoTest.check(mNano, PopulationNanoTest.DENOVO_REF, sam, PopulationNanoTest.FAMILY_PED, testName, new String[]{"--Xindel-trigger-fraction", "0.5"}, 1200, true);
  }
  public void testIndelTriggeringOnTwoSplitIndels() throws Exception {
    final String testName = "indel_triggering_on2_split";
    final String sam = FileHelper.resourceToString(RESOURCE_DIR + testName + ".sam");
    PopulationNanoTest.check(mNano, PopulationNanoTest.DENOVO_REF, sam, PopulationNanoTest.FAMILY_PED, testName, new String[]{}, 1000, true);
  }
  public void testIndelTriggeringOnTwoSame() throws Exception {
    final String testName = "indel_triggering_on2_same";
    final String sam = FileHelper.resourceToString(RESOURCE_DIR + testName + ".sam");
    PopulationNanoTest.check(mNano, PopulationNanoTest.DENOVO_REF, sam, PopulationNanoTest.FAMILY_PED, testName, new String[]{}, 1000, true);
  }
  public void testIndelTriggeringDeletionLength() throws Exception {
    final String testName = "indel_triggering_deletion";
    final String sam = FileHelper.resourceToString(RESOURCE_DIR + testName + ".sam");
    PopulationNanoTest.check(mNano, PopulationNanoTest.DENOVO_REF, sam, PopulationNanoTest.FAMILY_PED, testName, new String[]{"--Xcontrary-probability", "1"}, 1760, true);
  }
  public void testIndelTriggeringDeletionReverse() throws Exception {
    final String testName = "indel_triggering_deletion_reverse";
    final String sam = FileHelper.resourceToString(RESOURCE_DIR + testName + ".sam");
    PopulationNanoTest.check(mNano, PopulationNanoTest.DENOVO_REF, sam, PopulationNanoTest.FAMILY_PED, testName, new String[]{"--Xcontrary-probability", "1"}, 1760, true);
  }

  public void testStartOfReferenceN() throws Exception {
    // If you call something that needs previous base in the first template position don't output it
    final String testName = "startOfReferenceN";
    final String sam = FileHelper.resourceToString(RESOURCE_DIR + testName + ".sam");
    PopulationNanoTest.check(mNano, PopulationNanoTest.DENOVO_REF, sam, PopulationNanoTest.FAMILY_PED, testName, new String[]{}, 260, false);
  }

  public void testBedRegion() throws Exception {
    try (TestDirectory dir = new TestDirectory()) {
      final File bedFile = new File(dir, "testRegion.bed");
      FileUtils.copyResource(RESOURCE_DIR + "testBedRegion.bed", bedFile);
      final String testName = "testBedRegion";
      final String sam = FileHelper.resourceToString(RESOURCE_DIR + testName + ".sam");
      final String template = FileHelper.resourceToString(RESOURCE_DIR + "bed_region_template.fasta");
      PopulationNanoTest.check(mNano, template, sam, PopulationNanoTest.FAMILY_PED, testName, new String[]{"--Xsimple-repeat-extension=false", "--filter-bed", bedFile.getPath()}, 1000, false);
    }
  }
  public void testRegionPlusBedRegion() throws Exception {
    try (TestDirectory dir = new TestDirectory()) {
      final File bedFile = new File(dir, "testRegion.bed");
      FileUtils.copyResource(RESOURCE_DIR + "testBedRegion.bed", bedFile);
      final String testName = "testBedRegionPlusRegion";
      final String sam = FileHelper.resourceToString(RESOURCE_DIR + "testBedRegion.sam");
      final String template = FileHelper.resourceToString(RESOURCE_DIR + "bed_region_template.fasta");
      PopulationNanoTest.check(mNano, template, sam, PopulationNanoTest.FAMILY_PED, testName, new String[]{"--Xsimple-repeat-extension=false", "--filter-bed", bedFile.getPath(), "--region", "seq1:10-90"}, 600, false);
    }
  }
}
