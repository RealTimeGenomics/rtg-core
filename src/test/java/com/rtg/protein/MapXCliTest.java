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
package com.rtg.protein;

import static com.rtg.util.StringUtils.LS;
import static com.rtg.util.cli.CommonFlagCategories.REPORTING;
import static com.rtg.util.cli.CommonFlagCategories.SENSITIVITY_TUNING;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.MainResult;
import com.rtg.mode.DnaUtils;
import com.rtg.mode.SequenceMode;
import com.rtg.ngs.MapFlags;
import com.rtg.ngs.NgsFilterParams;
import com.rtg.ngs.NgsFilterParams.NgsFilterParamsBuilder;
import com.rtg.ngs.NgsParams;
import com.rtg.ngs.OutputFilter;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.IntegerOrPercentage;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.Utils;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Flag;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

/**
 * Test for corresponding class
 */
public class MapXCliTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new MapXCli();
  }

  public void testMapXFlags() {
    assertEquals("mapx", mCli.moduleName());
    checkHelp("directory for output",
        "SDF containing protein database to search",
        "guaranteed number of positions",
        "guaranteed minimum number of gaps",
        "scoring matrix",
        "maximum repeat frequency",
        "guaranteed minimum number of identical mismatches",
        "word size",
        "maximum alignment score at output",
        "maximum e-score at output",
        "maximum number of topn/topequals results output",
        "min-bit-score", "minimum bit score at output",
        "minimum percent identity at output",
        "output filter",
        "output unmapped reads",
        "do not gzip the output",
        "directory used for temporary files",
        "number of threads",
        "min-identity", "minimum percent identity at output"
        );
  }

  static final String TEMPLATE_DNA = "AAATGGCGCAAAAACAGAAAGTCGAAAAAAAATCAAAGAAATTATAACCACGACGCAGCAGACGCAG";

  static final String TEMPLATE_PROTEIN = TestUtils.dnaToProtein(TEMPLATE_DNA);

  static final String TEMPLATE_FASTA = ">templateName" + LS + TEMPLATE_PROTEIN + LS;

  private static final String[] READS_ONE_INDEL = {
    "TGGCG" + "ACG" + "CAAAAACAGAAAGTCGAAAAAAAATCAA",
    "AATGGCGCAAAAACAGAAAGTCATGGAAAAAAAATC",
    "ATGGCGCAAAAACACGCGAAAGTCGAAAAAAAATCA",
    "TGGCGCAAAAACAGAAAGTCGAAATTTAAAAATCAA",
    DnaUtils.reverseComplement("AAATGGCGCAAAAAACCCAGAAAGTCGAAAAAAAAT"),
    DnaUtils.reverseComplement("AATGGCGCAAAAACAGAAAGTTGACGAAAAAAAATC"),
    DnaUtils.reverseComplement("ATGGCGCAAAAACAGAAAGTCGGGGAAAAAAAATCA"),
    DnaUtils.reverseComplement("TGGCGCAAAAACAGAAAGTCGACACAAAATCACGAA"),
  "AAAAAAATCAAAGAAATTATAACCACGACAAAGCAG"};
  private static final String[] READS_MULTI_SUB = {
    "AAATGGCGCAAAAACAGACAGTCGAAAAAAAGTCAA",
    "AATGGCGCAAAGACAGAAAGTCGAAAAAACATCAAA",
    "ATGGCGCAAAAACAGAAAGTCGAAATTAAATCAAAG",
    "TGGCGCAAAAACAGAAAGTCGACAAAAAATCAAAGA",
    DnaUtils.reverseComplement("AAATGGCGCAAAAACAGAAAGTCGAAAAAAAATCAA"),
    DnaUtils.reverseComplement("AATGCAGCAAAAACAGAAAGTCGAAAAAAAATCAAA"),
    DnaUtils.reverseComplement("ATGGCCCAAAAACAGAAAGTCGAAAAAAAATCAAAG"),
    DnaUtils.reverseComplement("TGGCGGAATAACAGAAAGTCGAAAAAAAATCAAAGA"),
    "AAAAAAATCAAAGAAATTATAACCACGACACGGCAG"
  };
  static final String READS_FASTA_ONE_INDEL;
  static final String READS_FASTA_MULTI_SUB;

  static {
    StringBuilder sb = new StringBuilder();
    for (int i = 0; i < READS_ONE_INDEL.length; ++i) {
      sb.append(">testRead").append(i).append(LS).append(READS_ONE_INDEL[i]).append(LS);
    }
    READS_FASTA_ONE_INDEL = sb.toString();
    sb = new StringBuilder();
    for (int i = 0; i < READS_MULTI_SUB.length; ++i) {
      sb.append(">testRead").append(i).append(LS).append(READS_MULTI_SUB[i]).append(LS);
    }
    READS_FASTA_MULTI_SUB = sb.toString();
  }

  public void testFlags() throws IOException {
    final File dir = FileUtils.createTempDir("mapx", "test");
    try {
      final File output = new File(dir, "output");
      final File reads = new File(dir, "reads");
      final File template = new File(dir, "template");
      ReaderTestUtils.getReaderProtein(TEMPLATE_FASTA, template).close();
      ReaderTestUtils.getReaderDNA(READS_FASTA_ONE_INDEL, reads, null).close();
      checkHandleFlagsOut("-t", template.getPath(), "-i", reads.getPath(), "-a", "1", "-b", "1", "-o", output.getPath(), "-w", "4", "-T", "1");
      CFlags flags = getCFlags();
      assertEquals("blosum62", flags.getValue(MapXCli.MATRIX_FLAG));
      assertEquals(IntegerOrPercentage.valueOf("30%"), flags.getValue(MapXCli.MAX_ALIGNMENT_SCORE));
      assertEquals(IntegerOrPercentage.valueOf("95%"), flags.getValue(CommonFlags.REPEAT_FREQUENCY_FLAG));
      assertEquals(1, flags.getValue(MapXCli.GAP_LENGTH_FLAG));

      checkHandleFlagsOut("-t", template.getPath(), "-i", reads.getPath(), "-a", "1", "-b", "1", "-o", output.getPath(), "-w", "4", "-T", "1");
      checkHandleFlagsOut("-t", template.getPath(), "-i", reads.getPath(), "-a", "1", "-b", "1", "-o", output.getPath(), "-w", "4", "-T", "1", "--matrix", "blosum45");
      checkHandleFlagsOut("-t", template.getPath(), "-i", reads.getPath(), "-a", "1", "-b", "1", "-o", output.getPath(), "-w", "4", "-T", "1", "--matrix", "blosum62");
      checkHandleFlagsOut("-t", template.getPath(), "-i", reads.getPath(), "-a", "1", "-b", "1", "-o", output.getPath(), "-w", "4", "-T", "1", "--matrix", "blosum80", "-c", "10", "-e", "15");

      flags = getCFlags();
      assertEquals(10, flags.getValue(MapXCli.GAP_LENGTH_FLAG));
      assertEquals(IntegerOrPercentage.valueOf(15), flags.getValue(MapXCli.MAX_ALIGNMENT_SCORE));

      checkHandleFlagsOut("-t", template.getPath(), "-i", reads.getPath(), "-a", "1", "-b", "1", "-o", output.getPath(), "-w", "4", "-T", "1", "--matrix", "blosum80", "-c", "10", "-e", "15", "--no-unmapped");
      assertEquals(true, getCFlags().isSet(MapFlags.NO_UNMAPPED));

      checkHandleFlagsErr("-t", template.getPath(), "-i", reads.getPath(), "-a", "1", "-b", "1", "-o", output.getPath(), "-w", "4", "-T", "1", "--step");
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  public void testModuleName() {
    assertEquals("mapx", getCli().moduleName());
  }

  public void testFiltering()  throws InvalidParamsException, IOException  {
    Diagnostic.setLogStream();
    final File dir = FileUtils.createTempDir("mapx", "test");
    try {
      final File template = new File(dir, "template");
      final File reads = new File(dir, "reads");
      final File output = new File(dir, "output");
      ReaderTestUtils.getReaderProtein(TEMPLATE_FASTA, template).close();
      ReaderTestUtils.getReaderDNA(READS_FASTA_MULTI_SUB, reads, null).close();

      NgsParams p = createParams(new String[] {"-t", template.getPath(), "-i", reads.getPath(), "-a", "1", "-b", "0", "-o", output.getPath(), "-w", "9", "-Z", "--all-hits"});
      assertEquals("PROTEIN_ALL_HITS", p.outputParams().filter().outputFilter().name());

      p = createParams(new String[] {"-f", "topn", "-t", template.getPath(), "-i", reads.getPath(), "-a", "1", "-b", "0", "-o", output.getPath(), "-w", "9", "-E", "2"});
      assertEquals("PROTEIN_TOPN", p.outputParams().filter().outputFilter().name());
      assertEquals(2.0D, p.outputParams().filter().maxEScore());

      p = createParams(new String[] {"-f", "topequal", "-t", template.getPath(), "-i", reads.getPath(), "-a", "1", "-b", "0", "-o", output.getPath(), "-w", "9"});
      assertEquals("PROTEIN_TOPEQUAL", p.outputParams().filter().outputFilter().name());
      assertEquals(10.0D, p.outputParams().filter().maxEScore());

      p = createParams(new String[] {"-t", template.getPath(), "-i", reads.getPath(), "-a", "1", "-b", "0", "-o", output.getPath(), "-w", "9", "-B", "5.0"});
      assertEquals(Double.MAX_VALUE, p.outputParams().filter().maxEScore());
      assertEquals(5.0D, p.outputParams().filter().minBitScore());
      assertEquals(60, p.outputParams().filter().minIdentity());

      p = createParams(new String[] {"-t", template.getPath(), "-i", reads.getPath(), "-o", output.getPath(), "-w", "9", "-B", "5.0", "-P", "30"});
      assertEquals(Double.MAX_VALUE, p.outputParams().filter().maxEScore());
      assertEquals(30, p.outputParams().filter().minIdentity());
      assertEquals(1, p.maskParams().getSubstitutions());
      assertEquals(0, p.maskParams().getIndels());
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  public void testBuildFilter() {

    final CFlags flags = new CFlags("blah", TestUtils.getNullPrintStream(), TestUtils.getNullPrintStream());
    flags.registerOptional('e', MapXCli.MAX_ALIGNMENT_SCORE, IntegerOrPercentage.class, CommonFlags.INT, "maximum alignment score at output (as absolute value or percentage of read length)", IntegerOrPercentage.valueOf("10%")).setCategory(REPORTING);
    flags.registerOptional('n', MapFlags.MAX_TOP_RESULTS_FLAG, Integer.class, CommonFlags.INT, "maximum number of top equal results output per read", 5).setCategory(REPORTING);
    flags.registerOptional(MapFlags.XSCORE_INDEL, Integer.class, CommonFlags.INT, "set max score indel for topn threshold", MapFlags.MAX_SCORE).setCategory(REPORTING); //7 was used for illumina mappings
    flags.registerOptional(CommonFlags.EXCLUDE_FLAG, "any results with an E-value larger than the specified value are suppressed in the output").setCategory(CommonFlagCategories.UTILITY);
    flags.registerOptional(MapFlags.TOPN_RESULTS_FLAG, Integer.class, CommonFlags.INT, "set the number of results per read for topn. Allowed values are between 1 and 255", 5).setCategory(REPORTING);
    flags.registerOptional('P', MapXCli.MIN_IDENTITY_FLAG, Integer.class, CommonFlags.INT, "minimum percent identity at output", 60).setCategory(REPORTING);
    final Flag<String> filter = flags.registerOptional('f', CommonFlags.OUTPUT_FILTER_FLAG, String.class, "name", "output filter", "topn");
    filter.setParameterRange(MapXCli.FILTERS.keySet());
    flags.registerOptional('E', MapXCli.MAX_ESCORE_FLAG, Double.class, CommonFlags.FLOAT, "maximum e-score at output", 10.0).setCategory(REPORTING);
    flags.registerOptional('B', MapXCli.MIN_BITSCORE_FLAG, Double.class, CommonFlags.FLOAT, "minimum bit score at output").setCategory(REPORTING);
    flags.registerOptional(MapXCli.PRE_FILTER_ALGORITHM, Integer.class, CommonFlags.INT, "pre-filter algorithm", -3).setCategory(SENSITIVITY_TUNING);


    final NgsFilterParamsBuilder nfpb = new NgsFilterParams.NgsFilterParamsBuilder();

    flags.setFlags("-e", "3", "--Xexclude", "true");
    MapXCli.buildFilterParams(flags, nfpb);

    final NgsFilterParams nfp = nfpb.create();
    assertEquals(OutputFilter.PROTEIN_TOPN, nfp.outputFilter());
    assertEquals(10, nfp.maxTopResults());
    assertTrue(nfp.exclude());
    assertFalse(nfp.useids());
    assertEquals(new IntegerOrPercentage(3), nfp.matedMaxMismatches());
  }

  public void testParams() throws Exception  {
    final File dir = FileUtils.createTempDir("mapx", "test");
    try {
      final File template = new File(dir, "template");
      final File reads = new File(dir, "reads");
      final File output = new File(dir, "output");
      ReaderTestUtils.getReaderProtein(TEMPLATE_FASTA, template).close();
      ReaderTestUtils.getReaderDNA(READS_FASTA_MULTI_SUB, reads, null).close();
      final NgsParams p = createParams(new String[] {"-t", template.getPath(), "-i", reads.getPath(), "-a", "1", "-b", "0", "-o", output.getPath(), "-w", "9", "--Xdont-merge-alignment-result", "--matrix", "blosum45"});
      assertEquals(output, ((MapXCli) mCli).outputDirectory());
      assertEquals(template, p.searchParams().directory());
      //assertTrue(p.searchParams().readerParams().reader() instanceof CompressedMemorySequencesReader2);
      assertEquals(SequenceMode.PROTEIN, p.searchParams().mode());
      assertEquals(reads, p.buildFirstParams().directory());
      //assertTrue(p.buildFirstParams().readerParams().reader() instanceof CompressedMemorySequencesReader2);
      assertEquals(SequenceMode.TRANSLATED, p.buildFirstParams().mode());

      assertFalse(p.outputParams().mergeMatchResults());
      assertFalse(p.outputParams().mergeAlignmentResults());

      assertNotNull(p.proteinScoringMatrix());

      assertEquals(9, p.maskParams().getWordSize());

    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  private NgsParams createParams(final String[] params) throws InvalidParamsException, IOException {
    checkHandleFlagsOut(params);
    return ((MapXCli) mCli).makeParams();
  }

  static final String WARN_READS = ""
      + ">small" + LS
      + "acgtacgtacg" + LS
      + ">big" + LS
      + "acgtacgtacgt" + LS;

  public void testWarningAndError() throws IOException, InvalidParamsException {
    try (final TestDirectory dir = new TestDirectory("mapx")) {
      final MemoryPrintStream ps = new MemoryPrintStream();
      Diagnostic.setLogStream(ps.printStream());
      final File template = new File(dir, "template");
      final File warnReads = new File(dir, "wreads");
      final File output = new File(dir, "output");
      ReaderTestUtils.getReaderProtein(TEMPLATE_FASTA, template).close();
      ReaderTestUtils.getReaderDNA(WARN_READS, warnReads, null).close();
      createParams(new String[] {"-t", template.getPath(), "-i", warnReads.getPath(), "-a", "0", "-b", "0", "-o", output.getPath(), "-w", "1", "--min-read-length", "12"});
      assertTrue(ps.toString(), ps.toString().contains("The input set contains reads which are shorter than the minimum query length 12 which will be ignored"));
      try {
        createParams(new String[] {"-t", template.getPath(), "-i", warnReads.getPath(), "-a", "0", "-b", "0", "-o", output.getPath(), "-w", "1", "--min-read-length", "13"});
        fail();
      } catch (final NoTalkbackSlimException e) {
        assertEquals("All reads are shorter than the minimum query length 13", e.getMessage());
      }
    }
  }

  private static final String PROTEIN_TEMPLATE = ">template_sequence" + LS
      + "WSEITLAVMFAGKWK*HKESKIF*PLPTPSYHSVGPVAHSVCGEYPHCLFAVNSSGVVPTMHRVRLVLSAGNTSCGCGSA"
      + "*RHRSEYA*NL*TGPNGNRTASATEIDGSKNSGCGLTVSEKSEGRQIHTIPKALWPFVLDLQHITSNSTPPQAEMVWVFQ"
      + "AGYKFCPRL*LKYEPVCTRRHSQINVAHSVLPAKIYGGVDALLIIAMCSYGRVNILRSNKIRSSYSRALAMEGTMLIGLT"
      + "PYLSGRAEPNLSIRIWQ*SRVRQRSHFNKVGTQAPLTRVSRSGKELLTGLA*CNQI*DIYSEVSIVVQDSHRVHGRLHDS*IMGSAAPPANAT" + LS;

  private static final String R_1 = "GAGCAACTCCACTCCACCCCAAGCCGAGATGGTCTGGGTATTTCAAGCGGGGTACAAGTT";
  private static final String READ_1 = ">read9/0/0/simulatedSequence1/435/F/60." + LS
      + R_1 + LS;

  private static final String R_2 = "TAGGACGGCGTAGGGAGTGGCTAGAAAATCTTCGATTCCTTGTGCTATTTCCACTTACCA";
  private static final String READ_2 = ">read7/0/0/simulatedSequence1/33/R/60." + LS
      + R_2 + LS;

  private static final String R_3 = "TCTTCGCGGGCAAGACAGAGTGCGCAACGTTGATCTGTGAGTGTCTGCGTGTGCAAACGGGTTCATACTT";
  private static final String READ_3 = ">read1/0/0/simulatedSequence1/514/R/70." + LS
      + R_3 + LS;

  private static final String R_4 = "AGGACTGCATCGGCAACTGAGATCGACGGCAGCAAAAACTCAGGATGTGGTCTGACTGTGTCCGAAAAGTCCGAGGGACGTCAGATTC";
  private static final String READ_4 = ">read3/0/0/simulatedSequence1/295/F/88." + LS
      + R_4 + LS;
  private static final long LENGTH = R_1.length() + R_2.length() + R_3.length() + R_4.length();

  public void testVariableLength() throws Exception {
    checkVariableLength(new HashSet<>(Arrays.asList(0, 1, 2, 3)), LENGTH, "-Z", "-a", "2");
  }
  public void testMinReadLength() throws Exception {
    checkVariableLength(new HashSet<>(Arrays.asList(2, 3)), 158L, "-Z", "-a", "2", "--min-read-length", "70");
  }
  public void testStartRead() throws Exception {
    checkVariableLength(new HashSet<>(Arrays.asList(1, 2, 3)), 218L, "-Z", "-a", "2", "--start-read", "1");
  }
  public void testEndRead() throws Exception {
    checkVariableLength(new HashSet<>(Arrays.asList(0, 1, 2)), 190L, "-Z", "-a", "2", "--end-read", "3");
  }
  public void testReadRange() throws Exception {
    checkVariableLength(new HashSet<>(Arrays.asList(1, 2)), 130L, "-Z", "-a", "2", "--start-read", "1", "--end-read", "3");
  }
  public void testNoCompress() throws Exception {
    checkVariableLength(new HashSet<>(Arrays.asList(1, 2)), 130L, "-Z", "-a", "2", "--start-read", "1", "--end-read", "3", "--Xcompress-hashes", "false");
  }
  public void checkVariableLength(Set<Integer> expectedMappings, long usageMetric, String... args) throws IOException {
    try (final TestDirectory dir = new TestDirectory("mapx")) {
      final File output = new File(dir, "output");
      final File template = new File(dir, "template");
      ReaderTestUtils.getReaderProtein(PROTEIN_TEMPLATE, template).close();
      final File reads = ReaderTestUtils.getDNADir(READ_1 + READ_2 + READ_3 + READ_4, new File(dir, "reads"));

      final AbstractCli mapXCli = getCli();
      MainResult result = MainResult.run(mapXCli, Utils.append(args, "-i", reads.getPath(), "-t", template.getPath(), "-o", output.getPath()));
      assertEquals(result.err(), 0, result.rc());
      final String[] res = StringUtils.grep(FileUtils.fileToString(new File(output, "alignments.tsv")), "^[^#]").split("\n|\r\n");
      final HashSet<Integer> readsSet = new HashSet<>();
      for (final String s : res) {
        final String[] line = s.split("\t");
        assertEquals("template_sequence", line[0]);
        readsSet.add(Integer.parseInt(line[2]));
      }
      assertEquals(expectedMappings, readsSet);
      final String usageLog = mapXCli.usageLog();
      TestUtils.containsAll(usageLog, "[Usage beginning module=mapx runId=", ", Usage end module=mapx runId=", " metric=" + usageMetric + " success=true]");
    }
  }

  public void testLongReadMetaChunking() throws IOException {
    try (final TestDirectory dir = new TestDirectory("mapx")) {
      final File output = new File(dir, "output");
      final File template = new File(dir, "template");
      ReaderTestUtils.getReaderProtein(FileHelper.resourceToString("com/rtg/protein/resources/mcProt.fa"), template).close();
      final File reads = ReaderTestUtils.getDNADir(FileHelper.resourceToString("com/rtg/protein/resources/mcReads.fa"), new File(dir, "reads"));
      checkMainInitOk("-i", reads.getPath(), "-t", template.getPath(), "-o", output.getPath());
      mNano.check("mcResults", StringUtils.grep(FileHelper.gzFileToString(new File(output, "alignments.tsv.gz")), "^[^#]"));
    }
  }
}
