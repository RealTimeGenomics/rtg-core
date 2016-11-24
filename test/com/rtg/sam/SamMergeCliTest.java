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
package com.rtg.sam;


import static com.rtg.util.StringUtils.LS;
import static com.rtg.util.StringUtils.TAB;

import java.io.File;
import java.io.IOException;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.tabix.TabixIndexer;
import com.rtg.util.Environment;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

/**
 */
public class SamMergeCliTest extends AbstractCliTest {

  /** Record marked as pcr duplicates. */
  public static final String SAM_PCR_DUP = ""
      + "555" + TAB + "1171" + TAB + "g1" + TAB + "11" + TAB + "255" + TAB + "8M" + TAB + "=" + TAB + "3" + TAB + "0" + TAB + "TTCAGCTA" + TAB + "````````" + TAB + "AS:i:1" + LS;

  @Override
  protected AbstractCli getCli() {
    return new SamMergeCli();
  }

  public void testHelp() {
    checkHelp("rtg sammerge"
        , "coordinate-sorted SAM/BAM files."
        , "-I,", "--input-list-file=FILE", "file containing a list of SAM/BAM format files (1 per line) containing mapped reads"
        , "-o,", "--output=FILE", "output SAM/BAM file"
        , "FILE+", "SAM/BAM format files containing coordinate-sorted reads"
        , "\\0, --exclude-mated", "exclude all mated SAM records"
        , "\\0, --exclude-unmated", "exclude all unmated SAM records"
        , "-m,", "--max-as-mated=INT", "if set, ignore mated SAM records with an alignment score (AS attribute) that exceeds this value"
        , "-u,", "--max-as-unmated=INT", "if set, ignore unmated SAM records with an alignment score (AS attribute) that exceeds this value"
        , "-c,", "--max-hits=INT", "if set, ignore SAM records with an alignment count that exceeds this value"
        , "\\0, --region=STRING", "if set, only process SAM records within the specified range. The format is one of <sequence_name>, <sequence_name>:start-end or <sequence_name>:start+length"
        , "-Z,", "--no-gzip", "do not gzip the output"
        );
  }

  public void testValidator() throws IOException {
    try (final TestDirectory temp = new TestDirectory("validator")) {
      final File fake = new File(temp, "fake.txt.sam");
      FileUtils.stringToFile(fake.getPath() + LS, fake);
      TestUtils.containsAll(checkHandleFlagsErr("-o", "blahOutput").replaceAll("\\s+", " "), "No input files specified.");
      TestUtils.containsAll(checkHandleFlagsErr("-o", fake.getPath(), fake.getPath(), "-Z").replaceAll("\\s+", " "), "The file \"" + fake.getPath() + "\" already exists. Please remove it first or choose a different file");
      TestUtils.containsAll(checkHandleFlagsErr("-o", fake.getPath(), fake.getPath(), "-c", "0").replaceAll("\\s+", " "), "--max-hits must be at least 1");
      checkHandleFlagsOut("-o", fake.getPath(), "-I", fake.getPath());
    }
  }

  public void testSimpleMerge() throws Exception {
    try (final TestDirectory temp = new TestDirectory("simplemerge")) {
      final File inputFile = new File(temp, "input.txt");
      final File ffa = new File(temp, "alignmentsA.sam.gz");
      FileHelper.stringToGzFile(SharedSamConstants.SAM1, ffa);
      new TabixIndexer(ffa, TabixIndexer.indexFileName(ffa)).saveSamIndex();
      final File ffb = new File(temp, "alignmentsB.sam.gz");
      FileHelper.stringToGzFile(SharedSamConstants.SAM9, ffb);
      new TabixIndexer(ffb, TabixIndexer.indexFileName(ffb)).saveSamIndex();
      FileUtils.stringToFile(ffa.getPath() + LS + ffb.getPath() + LS, inputFile);

      final File outFile = new File(temp, "test.sam.gz");
      final String stdout = checkMainInitOk("-I", inputFile.getPath(), "-o", outFile.getPath(), "--region", "g1:1+5");
      final String expected = "SAM records read:    6" + LS
                            + "SAM records written: 6" + LS;
      assertEquals(expected, stdout);
      assertTrue(outFile.exists());
      final String output = FileHelper.gzFileToString(outFile);
      TestUtils.containsAll(output
          , "SO:coordinate"
          , "@PG\tID:rtg", "VN:" + Environment.getVersion(), "PN:rtg"
          , "g1\t3", "g1\t5");
      assertEquals(output, 12, output.split("\r\n|\n").length);
      mNano.check("simplestdout", TestUtils.stripSAMHeader(output), false);
      assertTrue(TabixIndexer.indexFileName(outFile).exists());
    }
  }

  public void testSimpleMergeStdOut() throws Exception {
    Diagnostic.setLogStream();
    try (final TestDirectory temp = new TestDirectory("simplemerge")) {
      final File inputFile = new File(temp, "input.txt");
      final File ffa = new File(temp, "alignmentsA.sam.gz");
      FileHelper.stringToGzFile(SharedSamConstants.SAM1, ffa);
      new TabixIndexer(ffa, TabixIndexer.indexFileName(ffa)).saveSamIndex();
      final File ffb = new File(temp, "alignmentsB.sam.gz");
      FileHelper.stringToGzFile(SharedSamConstants.SAM9, ffb);
      new TabixIndexer(ffb, TabixIndexer.indexFileName(ffb)).saveSamIndex();
      FileUtils.stringToFile(ffa.getPath() + LS + ffb.getPath() + LS, inputFile);

      final String stdout = checkMainInitOk("-I", inputFile.getPath(), "--region", "g1:1+5");
      TestUtils.containsAll(stdout
        , "SO:coordinate"
        , "@PG\tID:rtg", "VN:" + Environment.getVersion(), "PN:rtg"
        , "g1\t3", "g1\t5");
      mNano.check("simplestdout", TestUtils.stripSAMHeader(stdout), false);
      assertEquals(stdout, 12, stdout.split("\r\n|\n").length);

      final File header = new File(temp, "header.sam.gz");
      FileHelper.stringToGzFile(""
        + "@HD" + TAB + "VN:1.0" + TAB + "SO:coordinate" + LS
        + "@SQ" + TAB + "SN:g1" + TAB + "LN:20" + LS
        + "@SQ" + TAB + "SN:gempty" + TAB + "LN:0" + LS
        + "@RG" + TAB + "ID:RG1" + TAB + "SM:TEST-UPDATED" + TAB + "PL:ILLUMINA" + LS, header);
      final String stdout2 = checkMainInitOk("-I", inputFile.getPath(), "--region", "g1:1+5", "--Xalternate-sam-header", header.getPath());
      mNano.check("simplestdout-newheader", TestUtils.sanitizeSAMHeader(stdout2), false);
    }
  }

  private static final String EXPECTED_TABIX_SAM = ""
   + "52\t0\tchr6\t160560249\t255\t93=\t*\t0\t0\tACAGCAGCAGCAACAACAGTAACAACAGTAGAAATAAGTACAATAGCAATAACAACAGTAATAGCAACAGCAAAAACAATAGCAGCAGTAACA\t*\tAS:i:0\tNM:i:0\tIH:i:1" + LS
   + "94\t0\tchr6\t160560323\t255\t20=\t*\t0\t0\tAACAATAGCAGCAGTAACAA\t*\tAS:i:0\tNM:i:0\tIH:i:1" + LS
   + "78\t16\tchr6\t160560324\t255\t51=\t*\t0\t0\tACAATAGCAGCAGTAACAATAACAACAGCAATAGCAGCAACAACAGCAACA\t*\tAS:i:0\tNM:i:0\tIH:i:1" + LS
   + "87\t16\tchr6\t160560357\t255\t51=\t*\t0\t0\tGCAGCAACAACAGCAACAAGAAAAATGACAATAGCAGCAGCAACAACAGCA\t*\tAS:i:0\tNM:i:0\tIH:i:1" + LS;


  public void testTabixHandlingBug() throws IOException {
    try (final TestDirectory dir = new TestDirectory("test")) {
      final File failing = FileHelper.resourceToFile("com/rtg/sam/resources/failing.sam.gz", new File(dir, "failing.sam.gz"));
      FileHelper.resourceToFile("com/rtg/sam/resources/failing.sam.gz.tbi", new File(dir, "failing.sam.gz.tbi"));
      final String out = checkMainInitOk(failing.getPath(), "--region", "chr6:160560323+100");
      assertEquals(EXPECTED_TABIX_SAM, StringUtils.grep(out, "^[^@].*$"));

      final String outlegacy = checkMainInitOk(failing.getPath(), "--region", "chr6:160560323+100", "--legacy-cigars");
      mNano.check("checklegacy", TestUtils.stripSAMHeader(outlegacy), false);
    }
  }

  public void testUnmappedMerge() throws Exception {
    Diagnostic.setLogStream();
    try (final TestDirectory temp = new TestDirectory("unmapped_merge")) {
      CommandLine.setCommandArgs("test", "arguments");
      final File inputFile = new File(temp, "input.txt");
      final File ffa = new File(temp, "alignmentsA.sam.gz");
      FileHelper.stringToGzFile(SharedSamConstants.SAM1, ffa);
      new TabixIndexer(ffa, TabixIndexer.indexFileName(ffa)).saveSamIndex();
      final File ffb = new File(temp, "alignmentsB.sam.gz");
      FileHelper.stringToGzFile(SharedSamConstants.SAM9 + SharedSamConstants.SAM_UNMAPPED, ffb);
      new TabixIndexer(ffb, TabixIndexer.indexFileName(ffb)).saveSamIndex();
      final File ffc = new File(temp, "alignmentsC.sam.gz");
      FileHelper.stringToGzFile(SharedSamConstants.SAM9 + SAM_PCR_DUP, ffc);
      new TabixIndexer(ffc, TabixIndexer.indexFileName(ffc)).saveSamIndex();
      FileUtils.stringToFile(ffa.getPath() + LS
        + ffb.getPath() + LS
        + ffc.getPath() + LS, inputFile);

      final File outFile = new File(temp, "test.sam.gz");
      final String stdout = checkMainInitOk("-I", inputFile.getPath(), "-o", outFile.getPath());
      mNano.check("sfm-inc-stdout", stdout);
      assertTrue(outFile.exists());
      final String output = FileHelper.gzFileToString(outFile);
      TestUtils.containsAll(output, "SO:coordinate", "@PG\tID:rtg", "VN:" + Environment.getVersion(), "PN:rtg");
      final String outRecords = TestUtils.stripSAMHeader(output);
      mNano.check("sfm-inc-unmappeddups", outRecords, false);
      assertTrue(TabixIndexer.indexFileName(outFile).exists());
    }
  }

  public void testUnmappedStrip() throws Exception {
    Diagnostic.setLogStream();
    try (final TestDirectory temp = new TestDirectory("unmapped_strip")) {
      CommandLine.setCommandArgs("test", "arguments");
      final File inputFile = new File(temp, "input.txt");
      final File ffa = new File(temp, "alignmentsA.sam.gz");
      FileHelper.stringToGzFile(SharedSamConstants.SAM1, ffa);
      new TabixIndexer(ffa, TabixIndexer.indexFileName(ffa)).saveSamIndex();
      final File ffb = new File(temp, "alignmentsB.sam.gz");
      FileHelper.stringToGzFile(SharedSamConstants.SAM9 + SharedSamConstants.SAM_UNMAPPED, ffb);
      new TabixIndexer(ffb, TabixIndexer.indexFileName(ffb)).saveSamIndex();
      final File ffc = new File(temp, "alignmentsC.sam.gz");
      FileHelper.stringToGzFile(SharedSamConstants.SAM9 + SAM_PCR_DUP, ffc);
      new TabixIndexer(ffc, TabixIndexer.indexFileName(ffc)).saveSamIndex();
      FileUtils.stringToFile(ffa.getPath() + LS
        + ffb.getPath() + LS
        + ffc.getPath() + LS, inputFile);

      final File outFile = new File(temp, "test.sam.gz");
      final String stdout = checkMainInitOk("-I", inputFile.getPath(), "-o", outFile.getPath(), "--exclude-duplicates", "--exclude-unmapped");
      mNano.check("sfm-exc-stdout", stdout);
      assertTrue(outFile.exists());
      final String output = FileHelper.gzFileToString(outFile);
      TestUtils.containsAll(output, "SO:coordinate", "@PG\tID:rtg", "VN:" + Environment.getVersion(), "PN:rtg");
      final String outRecords = TestUtils.stripSAMHeader(output);
      mNano.check("sfm-exc-unmappeddups", outRecords, false);
      assertTrue(TabixIndexer.indexFileName(outFile).exists());
    }
  }
}
