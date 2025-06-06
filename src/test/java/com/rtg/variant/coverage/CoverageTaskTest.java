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
package com.rtg.variant.coverage;

import static com.rtg.sam.SharedSamConstants.IN3SAM3HEADER;
import static com.rtg.sam.SharedSamConstants.OK0;
import static com.rtg.sam.SharedSamConstants.OK1;
import static com.rtg.sam.SharedSamConstants.OK2;
import static com.rtg.sam.SharedSamConstants.OK3;
import static com.rtg.sam.SharedSamConstants.OK4;
import static com.rtg.sam.SharedSamConstants.OK5;
import static com.rtg.sam.SharedSamConstants.OUT_SAM;
import static com.rtg.sam.SharedSamConstants.REF_SEQS;
import static com.rtg.sam.SharedSamConstants.REF_SEQS_M;
import static com.rtg.sam.SharedSamConstants.SAM1;
import static com.rtg.sam.SharedSamConstants.SAM10;
import static com.rtg.sam.SharedSamConstants.SAM3;
import static com.rtg.sam.SharedSamConstants.SAM9;
import static com.rtg.sam.SharedSamConstants.SAMHEADER1;
import static com.rtg.sam.SharedSamConstants.SAM_BODY;
import static com.rtg.sam.SharedSamConstants.SAM_M;
import static com.rtg.util.StringUtils.LS;
import static com.rtg.util.StringUtils.TAB;

import java.io.File;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.MainResult;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.tabix.IndexUtils;
import com.rtg.tabix.TabixIndexer;
import com.rtg.util.Environment;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.Utils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

/**
 */
public class CoverageTaskTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new CoverageCli();
  }

  public void test1() throws Exception {
    final String[] args0 = {
        "-s", "0", "-Z"
    };
    checkBed(REF_SEQS, SAM1, args0, "1", "", 0, false, false);
  }

  public void test1smooth1() throws Exception {
    final String[] args0 = {
        "-s", "1", "-Z"
    };
    checkBed(REF_SEQS, SAM1, args0, "1-smooth1", "", 0, false, false);
  }

  public void test1RC() throws Exception {
    final String[] args0 = {
        "-s", "0"
    };
    checkBed(REF_SEQS, SAM3, args0, "1", "", 0, true, false);
  }

  private static final String SAM_CIGAR1 = SAMHEADER1
  + "0" + TAB + "0" + TAB + "g1" + TAB +  "3" + TAB + "255" + TAB + "2M2D6M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ATCGACTG" + TAB + "````````" + TAB + "AS:i:0" + LS
  + "0" + TAB + "0" + TAB + "g1" + TAB + "11" + TAB + "255" + TAB + "1S1=2I5N2="   + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ATCGAC" + TAB + "``````" + TAB + "AS:i:0" + LS
  + "";

  /** Test some fancy cigars. */
  public void testCigar1() throws Exception {
    final String[] args0 = {
        "-s", "0"
    };
    checkBed(REF_SEQS, SAM_CIGAR1, args0, "2", "", 0, true, false);
  }

  /** Test some fancy cigars. */
  public void testCigar1Smooth1() throws Exception {
    final String[] args0 = {
        "-s", "1"
    };
    checkBed(REF_SEQS, SAM_CIGAR1, args0, "2-smooth1", "", 0, true, false);
  }

  /** Test some fancy cigars. */
  public void testCigar1Smooth2() throws Exception {
    final String[] args0 = {
        "-s", "2"
    };
    checkBed(REF_SEQS, SAM_CIGAR1, args0, "2-smooth2", "", 0, true, false);
  }

  /** Test multiple sequences */
  public void test3() throws Exception {
    final String[] args0 = {"-Z", "-s", "0"};
    checkBed(REF_SEQS_M, SAM_M, args0, "3", "", 0, false, false);
  }

  /** Test multiple sequences */
  public void test3WithRegion() throws Exception {
    final String[] args0 = {"-Z", "-s", "0", "--region", "g2:6+2"};
    checkBed(REF_SEQS_M, SAM_M, args0, "3region", "", 0, false, false);
  }

  public void test3WithBadRegion() throws Exception {
    final String[] args0 = {"-Z", "-s", "0", "--region", "notanactualregion"};
    check(REF_SEQS_M, SAM_M, args0, null, new String[] {"notanactualregion", "not found"}, null, 1, false, null, false, false, null);
  }

  public void testIHFilter2() throws Exception {
    final String[] args0 = {
        "-c", "1", "-s", "0"
    };
    checkBed(REF_SEQS, SAM9, args0, "5", "", 0, true, false);
  }

  public void testIHValidation() throws Exception {
    final String[] args0 = {
        "-c", "1"
    };
    checkBed(REF_SEQS, SAM10, args0, null, "", "4 records skipped due to input filtering criteria", 0, true, null, false);
  }


  /**
   * Check that reads with <code>IH=n</code> contribute only <code>1/n</code>.
   */
  public void testIHfraction() throws Exception {
    final String[] args0 = {"-s", "0"};
    checkBed(REF_SEQS, SAM9, args0, "5IHfraction", "", 0, true, false);
  }

  public void test3Multifile() throws Exception {
    final String[] args0 = {"-s", "0"};
    // several files already globally sorted.
    checkMultifile(REF_SEQS, args0, 1,
        IN3SAM3HEADER + OK0 + OK1,
        IN3SAM3HEADER + OK2 + OK3,
        IN3SAM3HEADER + OK4 + OK5);
    // several files, each sorted, but needing a merge sort.
    checkMultifile(REF_SEQS, args0, 1,
        IN3SAM3HEADER + OK0 + OK4,
        IN3SAM3HEADER + OK1 + OK3,
        IN3SAM3HEADER + OK2 + OK5);
  }

  public void testIHFilter1() throws Exception {
    final String[] args0 = {
        "-c", "1"
    };
    checkBed(REF_SEQS, SAM9, args0, null, "", "4 records skipped due to input filtering criteria", 0, true, null, false);
  }

  private static final String REF_SEQS_STATS = ""
    + ">g1" + LS
    + "aa" + "atcg" + "actg" + "gtca" + "gcta" + "gg" + LS
    + ">gempty" + LS
    + LS
    + ">g2" + LS
    + "aa" + "atcn" + "nctg" + "gtca" + "gcta" + "gg" + "aaaaacccnnnngggttttt" + LS;

  /** Multiple Sequence SAM Records **/
  private static final String SAM_M2 = ""
    + SAMHEADER1
    + "@SQ" + TAB + "SN:g2" + TAB + "LN:40" + LS
    + SAM_BODY.substring(SAM_BODY.indexOf(LS) + LS.length())  // remove first line
    + SAM_BODY.replace("g1" + TAB + "11", "g2" + TAB + "10").replace("g1", "g2")
    ;

  public void testStats() throws Exception {
    final String[] args0 = {"-Z", "-s", "0"};
    check(REF_SEQS_STATS, SAM_M2, args0, "Stats", new String[]{""}, "", 0, false, null, true, false, "LevelsStatsOut");
  }

  public void testThreads() throws Exception {
    //This is to check when you have specified 1 thread, we do not use ThreaderSamMultifileIterator
    final String[] args0 = {"-Z", "-T", "1"};
    checkBed(REF_SEQS, SAM9, args0, null, "", "", 0, true, null, false);

    final String[] args1 = {"-Z", "-T", "2"};
    checkBed(REF_SEQS, SAM9, args1, null, "", "", 0, true, null, false);
  }

  private void checkBed(String refSeq, String sam, String[] args0, String expFile, String errorMsg, int errCode, boolean gzip, boolean checkOut) throws Exception {
    check(refSeq, sam, args0, expFile, new String[] {errorMsg}, "", errCode, gzip, null, checkOut, false, null);
  }

  private void checkBed(String refSeq, String sam, String[] args0, String expFile,
      String errorMsg, String logMsg, int errCode, boolean gzip, String[] erlines, boolean checkOut) throws Exception {
    check(refSeq, sam, args0, expFile, new String[]{errorMsg}, logMsg, errCode, gzip, erlines, checkOut, false, null);
  }

  private void check(String refSeq, String sam, String[] args0, String expFile, String[] errorMsgs, String logMsg,
                     int errCode, boolean gzip, String[] erlines, boolean checkOut, boolean tsv, String levelExp) throws Exception {
    Diagnostic.setLogStream();
    try (final TestDirectory tmpDir = new TestDirectory()) {
      final File output = new File(tmpDir, "output_dir");
      final File samFile = new File(tmpDir, OUT_SAM);
      FileUtils.stringToFile(sam, samFile);
      final File samFileGz = IndexUtils.ensureBlockCompressed(samFile);
      new TabixIndexer(samFileGz).saveSamIndex();
      final String outn = output.getPath();

      final File templ = ReaderTestUtils.getDNADir(refSeq, new File(tmpDir, "refSeqDir"));
      final String[] coverageArgs = new String[]{
        "-t", templ.getPath(),
        "-o", outn,
        samFileGz.getPath(),
        "--Xdisable-html-report"
      };
      final MainResult res = MainResult.run(getCli(), Utils.append(coverageArgs, args0));
      assertEquals(res.err(), errCode, res.rc());
      //System.err.println("out\n" + out.toString() + "end out");
      if (checkOut) {
        mNano.check("coveragetasktest-out.txt", res.out());
      }
      //System.err.println(errStr);
      TestUtils.containsAll(res.err(), errorMsgs);
      final String lout = FileUtils.fileToString(new File(outn, "coverage.log"));
      if (logMsg != null && !logMsg.equals("")) {
        assertTrue(lout, lout.contains(logMsg));
      }
      if ((errorMsgs.length == 0 || (errorMsgs.length == 1 && errorMsgs[0].equals(""))) && expFile != null) {
        final String fileName = (tsv ? CoverageParams.TSV_NAME : CoverageParams.BED_NAME) + (gzip ? ".gz" : "");
        final File file = new File(output, fileName);
        //System.err.println(output.getCanonicalPath() + ":" + fileName + ":" + file.length());
        final String result = gzip ? FileHelper.gzFileToString(file) : FileUtils.fileToString(file);
        if (gzip) {
          assertTrue(new File(output, fileName + ".tbi").exists());
          assertTrue(lout.contains("CoverageIndex"));
        }
        //System.out.println(CoverageParams.NAME + ".txt" + " contains:\n" + result);
        mNano.check("coveragetasktest" + expFile + ".txt", result.replaceAll("#.*" + LS, ""));
        TestUtils.containsAll(result
          , "#Version " + Environment.getVersion() + ", Coverage" + (tsv ? "" : " BED") + " output " + CoverageBedWriter.COVERAGE_OUTPUT_VERSION
          , "#RUN-ID\t"
        );
        if (tsv) {
          assertTrue(result, result.contains("#sequence\tposition\tunique-count\tambiguous-count\tscore" + LS));
        }
        final File summaryFile = new File(outn, "summary.txt");
        assertTrue(summaryFile.exists());
        final String summary = FileHelper.fileToString(summaryFile);
        TestUtils.containsAll(res.out(), summary);

        if (levelExp != null) {
          final File levelsFile = new File(outn, "levels.tsv");
          final String levels = FileHelper.fileToString(levelsFile);
          mNano.check("coveragetasktest" + levelExp + ".txt", levels.replaceAll("#.*" + LS, ""));
        }
      }
      if (erlines != null) {
        // check sequencer.errors output file
        final File sefile = new File(output, "sequencer.errors");
        final String result = FileUtils.fileToString(sefile);
        TestUtils.containsAll(result, erlines);
      }
    }
  }

  private static File saveSamTabixed(File name, String contents) throws Exception {
    FileUtils.stringToFile(contents, name);
    final File ret = IndexUtils.ensureBlockCompressed(name);
    new TabixIndexer(ret).saveSamIndex();
    return ret;
  }

  private void checkMultifile(final String refSeq, final String[] args0, final int expNum, final String maps1, final String maps2, final String maps3) throws Exception {
    Diagnostic.setLogStream();
    try (final TestDirectory tmpDir = new TestDirectory()) {
      final File output = new File(tmpDir, "output_dir");
      final File sam1 = saveSamTabixed(new File(tmpDir, OUT_SAM + "1"), maps1);
      final File sam2 = saveSamTabixed(new File(tmpDir, OUT_SAM + "2"), maps2);
      final File sam3 = saveSamTabixed(new File(tmpDir, OUT_SAM + "3"), maps3);
      final File templ = ReaderTestUtils.getDNADir(refSeq, new File(tmpDir, "refSeqDir"));

      final String[] coverageArgs = {
              "-t", templ.getPath(),
              "-o", output.getPath(),
              sam1.getPath(),
              sam2.getPath(),
              sam3.getPath(),
              "--Xdisable-html-report"
      };
      final MainResult res = checkMainInit(Utils.append(coverageArgs, args0));
      assertEquals(res.err(), 0, res.rc());
      final String result = FileHelper.gzFileToString(new File(output, CoverageParams.BED_NAME + ".gz"));
      mNano.check("coveragetasktest" + expNum + ".txt", result.replaceAll("#.*" + LS, ""));
      TestUtils.containsAll(result
        , "#Version " + Environment.getVersion() + ", Coverage BED output " + CoverageBedWriter.COVERAGE_OUTPUT_VERSION
        , "#RUN-ID\t"
      );
    }
  }

  /** Paired-end variant of SAM1 **/
  public static final String SAM_IH1_NH2 = ""
    + SAMHEADER1
    //left paired end
    + "1" + TAB +   "0" + TAB + "g1" + TAB +  "3" + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ATCGACTG" + TAB + "````````" + TAB + "AS:i:0" + TAB + "IH:i:1" + TAB + "NH:i:2" + LS
    + "2" + TAB +   "0" + TAB + "g1" + TAB +  "5" + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACTGTT" + TAB + "````````" + TAB + "AS:i:1" + TAB + "IH:i:1" + TAB + "NH:i:2" + LS
    ;
  public void testNHandIH() throws Exception {
    check(REF_SEQS, SAM_IH1_NH2, new String[] {"-s", "0"}, "ih1nh2", new String[] {}, "", 0, true, null, false, false, null);
  }

  /** Paired-end variant of SAM1 **/
  public static final String SAM_IH2_NH5 = ""
    + SAMHEADER1
    //left paired end
    + "1" + TAB +   "0" + TAB + "g1" + TAB +  "3" + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ATCGACTG" + TAB + "````````" + TAB + "AS:i:0" + TAB + "IH:i:2" + TAB + "NH:i:5" + LS
    + "3" + TAB +   "0" + TAB + "g1" + TAB +  "3" + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ATCGACTG" + TAB + "````````" + TAB + "AS:i:0" + TAB + "IH:i:2" + TAB + "NH:i:5" + LS
    + "2" + TAB +   "0" + TAB + "g1" + TAB +  "5" + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACTGTT" + TAB + "````````" + TAB + "AS:i:1" + TAB + "IH:i:2" + TAB + "NH:i:5" + LS
    + "4" + TAB +   "0" + TAB + "g1" + TAB +  "5" + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACTGTT" + TAB + "````````" + TAB + "AS:i:1" + TAB + "IH:i:2" + TAB + "NH:i:5" + LS
    ;
  public void testNH5andIH2() throws Exception {
    //this should give exactly same results as testNHandIH, as we doubled the coverage and Ih = 1 => ih = 2
    check(REF_SEQS, SAM_IH2_NH5, new String[] {"-s", "0"}, "ih1nh2", new String[] {}, "", 0, true, null, false, false, null);
  }

  private static final String SAM_HEADER_CLIP = ""
          + "@HD" + TAB + "VN:1.0" + TAB + "SO:coordinate" + LS
          + "@SQ" + TAB + "SN:t" + TAB + "LN:84" + LS
      ;
  private static final String SAM_CLIP = SAM_HEADER_CLIP
          + "0" + TAB + "0" + TAB + "t" + TAB + "1" + TAB + "255" + TAB + "10S43M2S" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "AAAAAAAAAATCGCTAGGTTCGACTTGGTTAACAACAACGCCTGGGGCTTTTTGG" + TAB + "*" + LS
          + "1" + TAB + "0" + TAB + "t" + TAB + "42" + TAB + "255" + TAB + "10S40M2S" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "AAAAAAAAAATTATTCTGGAAAGCAATGCCAGGCAGGGGCAGGTGGCCACGG" + TAB + "*" + LS;


  public void testSoftClipping() throws Exception {
    try (final TestDirectory tmpDir = new TestDirectory()) {
      final File template = ReaderTestUtils.getDNADir(">t\nTCGCTAGGTTCGACTTGGTTAACAACAACGCCTGGGGCTTTTTATTCTGGAAAGCAATGCCAGGCAGGGGCAGGTGGCCACTTT\n", new File(tmpDir, "template"));
      final File samFile = IndexUtils.ensureBlockCompressed(FileUtils.stringToFile(SAM_CLIP, new File(tmpDir, "sam.sam")));
      new TabixIndexer(samFile).saveSamIndex();
      final File output = new File(tmpDir, "output");
      checkMainInitOk("-t", template.getPath(), "-o", output.getPath(), "-s", "0", samFile.getPath(), "--Xdisable-html-report");
      final String was = StringUtils.grepMinusV(FileHelper.gzFileToString(new File(output, CoverageParams.BED_NAME + FileUtils.GZ_SUFFIX)), "^#");
      mNano.check("coveragetasktest_clipping.bed", was.replaceAll(LS, "\n")); //Bed does not care of EOL chars
    }
  }

  public void testSoftClippingPerBase() throws Exception {
    try (final TestDirectory tmpDir = new TestDirectory()) {
      final File template = ReaderTestUtils.getDNADir(">t\nTCGCTAGGTTCGACTTGGTTAACAACAACGCCTGGGGCTTTTTATTCTGGAAAGCAATGCCAGGCAGGGGCAGGTGGCCACTTT\n", new File(tmpDir, "template"));
      final File samFile = IndexUtils.ensureBlockCompressed(FileUtils.stringToFile(SAM_CLIP, new File(tmpDir, "sam.sam")));
      new TabixIndexer(samFile).saveSamIndex();
      final File output = new File(tmpDir, "output");
      checkMainInitOk("-t", template.getPath(), "-o", output.getPath(), "--per-base", samFile.getPath(), "--Xdisable-html-report");
      final String was = StringUtils.grepMinusV(FileHelper.gzFileToString(new File(output, CoverageParams.TSV_NAME + FileUtils.GZ_SUFFIX)), "^#");
      mNano.check("coveragetasktest_clipping.tsv", was.replaceAll(LS, "\n")); //Bed does not care of EOL chars
    }
  }

  static final String BED_REGIONS = "simulatedSequence1\t40\t47\n"
                                            + "simulatedSequence1\t55\t65\tblah\n"
                                            + "simulatedSequence1\t60\t75\tfeh\n"
                                            + "simulatedSequence2\t17\t43\n"
                                            + "simulatedSequence2\t55\t59\tfeh\n"
                                            + "simulatedSequence2\t63\t94\n";

  public void testBedRegions() throws Exception {
    try (final TestDirectory tmpDir = new TestDirectory()) {
      final File samFile = new File(tmpDir, "sam.sam.gz");
      IndexUtils.ensureBlockCompressed(FileHelper.resourceToFile("com/rtg/variant/resources/coverage_mated.sam.gz", samFile));
      new TabixIndexer(samFile).saveSamIndex();
      final File bedRegionsFile = FileHelper.stringToGzFile(BED_REGIONS, new File(tmpDir, "bedRegions.bed.gz"));
      final File output = new File(tmpDir, "output");
      final String tmpl = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAANAAAAAAAAAANAAAAAAAAAAAAAAAAAAAAAAAAAA";
      final File template = ReaderTestUtils.getDNADir(">simulatedSequence1\n" + tmpl + "\n>simulatedSequence2\n" + tmpl + "\n", new File(tmpDir, "template"));
      checkMainInitWarn("-o", output.getPath(), "-s", "0", samFile.getPath(),
        "--bed-regions", bedRegionsFile.getPath(), "-t", template.getPath(), "--Xdisable-html-report");
      final String is = StringUtils.grepMinusV(FileHelper.gzFileToString(new File(output, CoverageParams.BED_NAME + FileUtils.GZ_SUFFIX)), "^#");
      mNano.check("covBedRegion.bed", is);
    }
  }

  public void testBedRegionsPerRegion() throws Exception {
    try (final TestDirectory tmpDir = new TestDirectory()) {
      final File samFile = new File(tmpDir, "sam.sam.gz");
      IndexUtils.ensureBlockCompressed(FileHelper.resourceToFile("com/rtg/variant/resources/coverage_mated.sam.gz", samFile));
      new TabixIndexer(samFile).saveSamIndex();
      final File bedRegionsFile = FileHelper.stringToGzFile(BED_REGIONS, new File(tmpDir, "bedRegions.bed.gz"));
      final File output = new File(tmpDir, "output");
      final String tmpl = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAANAAAAAAAAAANAAAAAAAAAAAAAAAAAAAAAAAAAA";
      final File template = ReaderTestUtils.getDNADir(">simulatedSequence1\n" + tmpl + "\n>simulatedSequence2\n" + tmpl + "\n", new File(tmpDir, "template"));
      checkMainInitWarn("-o", output.getPath(), samFile.getPath(),
        "--bed-regions", bedRegionsFile.getPath(), "--per-region", "-t", template.getPath(), "--Xdisable-html-report");
      final String is = StringUtils.grepMinusV(FileHelper.gzFileToString(new File(output, CoverageParams.BED_NAME + FileUtils.GZ_SUFFIX)), "^#");
      mNano.check("covBedRegionPerRegion.bed", is);
    }
  }

  public void testBedRegionsTsv() throws Exception {
    try (final TestDirectory tmpDir = new TestDirectory()) {
      final File samFile = new File(tmpDir, "sam.sam.gz");
      IndexUtils.ensureBlockCompressed(FileHelper.resourceToFile("com/rtg/variant/resources/coverage_mated.sam.gz", samFile));
      new TabixIndexer(samFile).saveSamIndex();
      final File bedRegionsFile = FileHelper.stringToGzFile(BED_REGIONS, new File(tmpDir, "bedRegions.bed.gz"));
      final File output = new File(tmpDir, "output");
      final String tmpl = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAANAAAAAAAAAANAAAAAAAAAAAAAAAAAAAAAAAAAA";
      final File template = ReaderTestUtils.getDNADir(">simulatedSequence1\n" + tmpl + "\n>simulatedSequence2\n" + tmpl + "\n", new File(tmpDir, "template"));
      checkMainInitWarn("-o", output.getPath(), samFile.getPath(), "--per-base",
        "--bed-regions", bedRegionsFile.getPath(), "-t", template.getPath(), "--Xdisable-html-report");
      final String is = StringUtils.grepMinusV(FileHelper.gzFileToString(new File(output, CoverageParams.TSV_NAME + FileUtils.GZ_SUFFIX)), "^#");
      mNano.check("covBedRegion.tsv", is);
    }
  }

  public void testBedRegionsTsvNoTemplate() throws Exception {
    try (final TestDirectory tmpDir = new TestDirectory()) {
      final File samFile = new File(tmpDir, "sam.sam.gz");
      IndexUtils.ensureBlockCompressed(FileHelper.resourceToFile("com/rtg/variant/resources/coverage_mated.sam.gz", samFile));
      new TabixIndexer(samFile).saveSamIndex();
      final File bedRegionsFile = FileHelper.stringToGzFile(BED_REGIONS, new File(tmpDir, "bedRegions.bed.gz"));
      final File output = new File(tmpDir, "output");
      checkMainInitWarn("-o", output.getPath(), samFile.getPath(), "--per-base",
        "--bed-regions", bedRegionsFile.getPath(), "--Xdisable-html-report");
      final File summary = new File(output, CommonFlags.SUMMARY_FILE);
      mNano.check("covBedRegionNoTemplate.txt", FileUtils.fileToString(summary));
    }
  }
}
