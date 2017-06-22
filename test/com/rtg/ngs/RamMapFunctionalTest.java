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
package com.rtg.ngs;

import static com.rtg.util.StringUtils.LS;

import java.io.File;
import java.io.IOException;

import com.rtg.alignment.EditDistanceFactory;
import com.rtg.launcher.AbstractNanoTest;
import com.rtg.launcher.MainResult;
import com.rtg.launcher.globals.CoreGlobalFlags;
import com.rtg.reader.FormatCli;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SdfId;
import com.rtg.util.PropertiesUtils;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.Utils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;


/**
 */
public class RamMapFunctionalTest extends AbstractNanoTest {

  private static final String TEMPLATE = ">template" + LS
  //1234567890123456789012345678901234567890
  + "tacgtcatgactgctgcanactgcatgcatgactg"
  + "actgactgcatatgcattcatactcatgatgcatgctgca"
  + "tgatctgcatactatcattttcatgcagtactgcatgcat"
  + "gcatgcatgctactgtacgatcgatcgatcgatcg" + LS
          + ">blanko" + LS + LS
  + ">t2" + LS
          + "acgtacgtacgtacgt";

  private static final String READ_LEFT = ">left" + LS
  + "tacgtcatgactgctgcatactgcatgcatgactg" + LS
  + ">leftunmapped" + LS
  + "tttttttttgggggggggggaaaaaaaaccccccc" + LS;
  private static final String READ_RIGHT = ">right" + LS
  + "gcatgcatgctactgtacgatcgatcgatcgatcg" + LS
   + ">rightunmapped" + LS
   + "tttttttttgggggggggggaaaaaaaaccccccc" + LS;

  public void testE2eNoMerge() throws Exception {
    end2End(false, 0, 34);
  }

  public void testE2eMerge() throws Exception {
    end2End(true, 0, 34);
  }

  public void testE2eNoMergeUnknownsPenalty() throws Exception {
    end2End(false, 2, 30);
  }

  public void end2End(boolean merge, int unknownsPenalty, int mapQ) throws IOException {
    try (final TestDirectory outer = new TestDirectory("map-end2end")) {
      final File template = new File(outer, "template");
      final File reads = new File(outer, "Reads");
      final File left = new File(reads, "left");
      final File right = new File(reads, "right");
      final File out = new File(outer, "out");

      final File header = new File(outer, "header");
      FileUtils.stringToFile("@RG\tID:L23\tSM:NA123", header);
      ReaderTestUtils.getReaderDNA(TEMPLATE, template, null).close();
      ReaderTestUtils.getReaderDNA(READ_LEFT, left, new SdfId(123L)).close();
      ReaderTestUtils.getReaderDNA(READ_RIGHT, right, new SdfId(123L)).close();

      String[] args = {
        "-i", reads.getPath(), "-t", template.getPath(), "-o", out.getPath(),
        "-e", "50", "-m", "0", "-w", "5",
        "--sam-rg", header.getPath(),
        "--" + MapFlags.GAP_OPEN_PENALTY_FLAG, "1",
        "--" + MapFlags.GAP_EXTEND_PENALTY_FLAG, "1",
        "--" + MapFlags.MISMATCH_PENALTY_FLAG, "1",
        "--" + MapFlags.UNKNOWNS_PENALTY_FLAG, "" + unknownsPenalty,
        "--" + MapFlags.SAM_FLAG,
        "--" + MapFlags.ALIGNER_BAND_WIDTH_FACTOR_FLAG, "1.0",
        "--aligner-mode", "general",
        "--XX" + CoreGlobalFlags.EDIT_DIST_MISMATCH_SOFT_CLIP, "0",
      };
      if (!merge) {
        args = Utils.append(args, "--" + MapFlags.DONT_UNIFY_FLAG);
      }

      final MainResult r = MainResult.run(new MapCli(), args);
      assertEquals(r.err(), 0, r.rc());

      final String mated = FileHelper.gzFileToString(new File(out, merge ? "alignments.sam.gz" : "mated.sam.gz"));
      TestUtils.containsAll(mated, "@HD",
        "@RG\tID:L23\tSM:NA123",
        "@PG\tID:rtg",
        "@CO\tTEMPLATE-SDF-ID:",
        "@CO\tREAD-SDF-ID:",
        "@SQ\tSN:template\tLN:150",
        "0\t67\ttemplate\t1\t" + mapQ + "\t18=1X16=",
        "0\t131\ttemplate\t116\t" + mapQ + "\t35=",
        "AS:i:" + unknownsPenalty + "\tNM:i:1\tXA:i:" + unknownsPenalty + "\tRG:Z:L23\tIH:i:1",
        "AS:i:0\tNM:i:0\tXA:i:" + unknownsPenalty + "\tRG:Z:L23\tIH:i:1");
      //System.err.println(mated);
      if (!merge) {
        final String unmated = FileHelper.gzFileToString(new File(out, "unmated.sam.gz"));
        TestUtils.containsAll(unmated, "@HD",
          "@RG\tID:L23\tSM:NA123",
          "@PG\tID:rtg",
          "@CO\tTEMPLATE-SDF-ID:",
          "@CO\tREAD-SDF-ID:",
          "@SQ\tSN:template\tLN:150"
        );
        //System.err.println(unmated);
      }
    }
  }

  static final String PARTIAL_TEMPLATE = ">template" + LS
  + "acgttgtgacgtacgtatgtgtagtagtgcgtcttactagtacccagtctgcgtcagactagactacgtctagacccatgacctgagctttgactgtgtg";


  static final String READS_LEFT =
    ">left1" + LS
    + "gttgt" + LS
    + ">left2" + LS
    + "cttac" + LS
    + ">left3" + LS
    + "ctaga" + LS;

  static final String READS_RIGHT =
    ">right1" + LS
    + "gtgcg" + LS
    + ">right2" + LS
    + "tctgc" + LS
    + ">right3" + LS
    + "gcttt" + LS;


  public void testEnd2EndPartial() throws IOException {
    try (final TestDirectory outer = new TestDirectory("map-end2end")) {
      final File template = new File(outer, "template");
      final File reads = new File(outer, "Reads");
      final File left = new File(reads, "left");
      final File right = new File(reads, "right");
      final File out = new File(outer, "out");

      final File header = new File(outer, "header");
      FileUtils.stringToFile("@RG\tID:L23\tSM:NA123", header);
      ReaderTestUtils.getReaderDNA(PARTIAL_TEMPLATE, template, null).close();
      ReaderTestUtils.getReaderDNA(READS_LEFT, left, new SdfId(123L)).close();
      ReaderTestUtils.getReaderDNA(READS_RIGHT, right, new SdfId(123L)).close();

      final MainResult r = MainResult.run(new MapCli(), "-i", reads.getPath(),
        "-t", template.getPath(),
        "-o", out.getPath(),
        "-E", "100%", "-m", "0",
        "--sam-rg", header.getPath(),
        "--start-read", "1",
        "--end-read", "2",
        "--" + MapFlags.SAM_FLAG,
        "--" + MapFlags.DONT_UNIFY_FLAG);
      assertEquals(r.err(), 0, r.rc());

      final String mated = FileHelper.gzFileToString(new File(out, "mated.sam.gz"));
      assertFalse(mated.contains("\n0"));
      assertFalse(mated.contains("\n2"));
      assertTrue(mated.contains("\n1"));
      TestUtils.containsAll(mated, "@HD",
          "@RG\tID:L23\tSM:NA123",
          "@PG\tID:rtg",
          "@CO\tTEMPLATE-SDF-ID:",
          "@CO\tREAD-SDF-ID:",
          "@SQ\tSN:template\tLN:100",
      "AS:i:0\tNM:i:0\tXA:i:0\tRG:Z:L23\tIH:i:1");
      //System.err.println(mated);
      final String unmated = FileHelper.gzFileToString(new File(out, "unmated.sam.gz"));
      TestUtils.containsAll(unmated, "@HD",
          "@RG\tID:L23\tSM:NA123",
          "@PG\tID:rtg",
          "@CO\tTEMPLATE-SDF-ID:",
          "@CO\tREAD-SDF-ID:",
          "@SQ\tSN:template\tLN:100"
      );
      //System.err.println(unmated);
    }
  }

  public void testEnd2EndAllHits() throws IOException {
    try (final TestDirectory outer = new TestDirectory("map-end2end")) {
      final File template = new File(outer, "template");
      final File reads = new File(outer, "Reads");
      final File left = new File(reads, "left");
      final File right = new File(reads, "right");
      final File out = new File(outer, "out");

      final File header = new File(outer, "header");
      FileUtils.stringToFile("@RG\tID:L23\tSM:NA123\tPL:IONTORRENT", header);
      ReaderTestUtils.getReaderDNA(PARTIAL_TEMPLATE, template, null).close();
      ReaderTestUtils.getReaderDNA(READS_LEFT, left, new SdfId(123L)).close();
      ReaderTestUtils.getReaderDNA(READS_RIGHT, right, new SdfId(123L)).close();

      final MainResult r = MainResult.run(new MapCli(), "-i", reads.getPath(),
        "-t", template.getPath(),
        "-o", out.getPath(),
        "--sam-rg", header.getPath(),
        "--start-read", "1",
        "--end-read", "2",
        "--all-hits",
        "--" + MapFlags.SAM_FLAG,
        "--" + MapFlags.DONT_UNIFY_FLAG);
      assertEquals(r.err(), 0, r.rc());
      final String mated = FileHelper.gzFileToString(new File(out, "alignments.sam.gz"));
      assertFalse(mated.contains("\n0"));
      assertFalse(mated.contains("\n2"));
      assertTrue(mated.contains("\n1"));
      TestUtils.containsAll(mated, "@HD",
        "@RG\tID:L23\tSM:NA123\tPL:IONTORRENT",
        "@PG\tID:rtg",
        "@CO\tTEMPLATE-SDF-ID:",
        "@CO\tREAD-SDF-ID:",
        "@SQ\tSN:template\tLN:100"
      );
      assertEquals("all-hits is set, quality calibration output is disabled." + LS + "all-hits is set, svprep output is disabled." + LS, r.err());
    }
  }

  public void testBug1351() throws Exception {
    try (final TestDirectory outer = new TestDirectory("map-end2end")) {
      final File template = new File(outer, "template");
      final File reads = new File(outer, "reads");
      final File out = new File(outer, "out");

      final String read1 = "TATATCGCAAGTTAATATCTGTTTGGAGGGTATCGGGTCAAGTTGGACGCGTTGCAGCGTGCCAGCCAGCATAGCCTATTAGAGATTGTCTGAGTTAACT";
      final String read2 = "ATCTTGTATATCGCAAGTTAATATCTGTTTGGAGGGTATCGGGTCAAGTTGGACGCGTTGCAGCGTGCCAGCCAGCATAGCCTATTAGAGATTGTCTGAG";
      final String templateStr = ">1" + LS + "nnnnnnnnnn" + read1 + "nnnnnnnnnn" + LS + ">2" + LS + "nnnnnnnnnn" + read2 + "nnnnnnnnnn" + LS;

      final String readStr = ">1" + LS + read1 + LS + ">2" + LS + read2 + LS;

      ReaderTestUtils.getReaderDNA(templateStr, template, null).close();
      ReaderTestUtils.getReaderDNA(readStr, reads, null).close();

      assertTrue(reads.canRead());

      final MainResult r = MainResult.run(new MapCli(), "-i", reads.getPath(),
        "-t", template.getPath(),
        "-o", out.getPath(),
        "--all-hits",
        "--soft-clip-distance", "0",
        "--" + MapFlags.SAM_FLAG,
        "--" + MapFlags.DONT_UNIFY_FLAG,
        "--" + MapFlags.UNKNOWNS_PENALTY_FLAG, "0",
        "--" + MapFlags.ALIGNER_MODE_FLAG, "general",
        "--XX" + CoreGlobalFlags.EDIT_DIST_MISMATCH_SOFT_CLIP, "2");
      assertEquals(r.err(), 0, r.rc());

      final String sam = FileHelper.gzFileToString(new File(out, "alignments.sam.gz"));

      final String[] split = StringUtils.grep(sam, "^[^@]").split(LS);
      assertEquals(sam, 4, split.length);
      assertTrue(sam, split[0].startsWith("0\t0\t1\t11"));
      assertTrue(sam, split[1].startsWith("1\t0\t1\t11"));
      assertTrue(sam, split[2].startsWith("1\t0\t2\t11"));
      assertTrue(sam, split[3].startsWith("0\t0\t2\t17"));
    }
  }

  public void testIonTorrentPair() throws Exception {
    try (final TestDirectory outer = new TestDirectory("map-end2end")) {
      final String templateStr = ">t" + LS + "GGTCCGTCACCTCCTGGTATTCCAGTTTGATGTTATCGGCAATACACAGCGCATTCATCAGTGGCTGGACGGTTTTTTCCTGCATATCCACGAGGTAGTAGTAACTGCCATCCTCGCAGCCATCGGGCGACTGGCGGGTACGTAAAACCTGCAAGGTGGGACCGGATAAGTAGTCAACCTGTGACCATTCTTCGCTGACATCATCCTGGTGGCTGACCACATCACTGAAACGCGGTGGGGTGAGGTCCTTAAATTTACTGATGGCCTTCAGGTCATCACTTCTGTCATCGCATGCGATGA" + LS;
      final String read1 = ">left" + LS + "GCCAGTCGCCCGATGGCTGCGAGGATGGCAGTTACTACTACCTCGTGGATATGCAGGAAAAAACCGTCCAGCCACTGATGAATGCGCTGTGTATTGCCGATAACATCAAGCTGAGGA" + LS;
      //    final String read1 = ">left" + StringUtils.LS + "TCCTCAGCTTGATGTTATCGGCAATACACAGCGCATTCATCAGTGGCTGGACGGTTTTTTCCTGCATATCCACGAGGTAGTAGTAACTGCCATCCTCGCAGCCATCGGGCGACTGGC" + StringUtils.LS;
      final String read2 = ">righ" + LS + "TTGATGTTATCGGCAATACACAGCGCATTCATCAGTGGCTGGACGGTTTTTTCCTGCATATCCACGAGGTAGTAGTAACTGCCATCCTCGCAGCCATCGGGCGACTGGCCTGAGTCGGAGACAC" + LS;

      final File template = new File(outer, "template");
      final File reads = new File(outer, "Reads");
      final File left = new File(reads, "left");
      final File right = new File(reads, "right");
      final File out = new File(outer, "out");

      ReaderTestUtils.getReaderDNA(templateStr, template, null).close();
      ReaderTestUtils.getReaderDNA(read1, left, new SdfId(123L)).close();
      ReaderTestUtils.getReaderDNA(read2, right, new SdfId(123L)).close();

      final MainResult r = MainResult.run(new MapCli(), "-i", reads.getPath(),
        "-t", template.getPath(),
        "-o", out.getPath(),
        //          "--sam-rg", header.getPath(),
        //          "--start-read", "1",
        //          "--end-read", "2"
        "--" + MapFlags.GAP_OPEN_PENALTY_FLAG, "1",
        "--" + MapFlags.GAP_EXTEND_PENALTY_FLAG, "1",
        "--" + MapFlags.MISMATCH_PENALTY_FLAG, "1",
        "--" + MapFlags.SOFT_CLIP_DISTANCE_FLAG, "0",
        "--" + MapFlags.SAM_FLAG,
        "--" + MapFlags.DONT_UNIFY_FLAG,
        "--" + MapFlags.ALIGNER_MODE_FLAG, "general",
        "--XX" + CoreGlobalFlags.EDIT_DIST_MISMATCH_SOFT_CLIP, "0");
      assertEquals(r.err(), 0, r.rc());

      final String mated = FileHelper.gzFileToString(new File(out, "mated.sam.gz"));
      TestUtils.containsAll(mated,
        "0\t83\tt\t20\t55\t1=1X1=1I3=1X109=\t=\t27\t129\tTCCTCAGCTTGATGTTATCGGCAATACACAGCGCATTCATCAGTGGCTGGACGGTTTTTTCCTGCATATCCACGAGGTAGTAGTAACTGCCATCCTCGCAGCCATCGGGCGACTGGC\t",
        "0\t163\tt\t27\t55\t109=2I1=1X2=2X1=2X1=1X2=\t=\t20\t-129\tTTGATGTTATCGGCAATACACAGCGCATTCATCAGTGGCTGGACGGTTTTTTCCTGCATATCCACGAGGTAGTAGTAACTGCCATCCTCGCAGCCATCGGGCGACTGGCCTGAGTCGGAGACAC\t"
      );
    }
  }

//  public void testRegression() {
//    final File dir = FileUtils.createTempDir("rammap", "functional");
//    final MemoryPrintStream mps = new MemoryPrintStream();
//    try {
//      final File out = new File(dir, "out");
//
//      final RamMap map = new RamMap();
//      final String[] args = {
//          "-i", "Y:\\users\\mehul\\workspace\\rtg\\internal\\com\\rtg\\regression\\data\\long_pe_rammap\\bug997\\r",
//          "-t", "Y:\\users\\mehul\\workspace\\rtg\\internal\\com\\rtg\\regression\\data\\long_pe_rammap\\bug997\\g",
//          "-o", out.getPath(),
//          "-a", "0",
//          "-b", "0",
//          "-m", "2000",
//          "-M", "2500",
////          "--sam-rg", header.getPath(),
////          "--start-read", "1",
////          "--end-read", "2"
//      };
//
//      final int code = map.mainInit(args, mps.outputStream(), mps.printStream());
//      assertEquals(mps.toString(), 0, code);
//
//      /*final String expected = FileUtils.fileToString("Y:\\users\\mehul\\workspace\\rtg\\internal\\com\\rtg\\regression\\data\\long_pe_rammap\\bug997\\expect_out\\unmapped.sam");
//      final String unmapped = FileHelper.gzFileToString(new File(out, "unmapped.sam.gz"));
//      assertEquals(expected.replaceAll("^@.*$", ""), unmapped.replaceAll("^@.*$", ""));
//*/
//      final String expected_unmated = FileUtils.fileToString("Y:\\users\\mehul\\workspace\\rtg\\internal\\com\\rtg\\regression\\data\\long_pe_rammap\\bug997\\expect_out\\unmated.sam");
//      final String unmated = FileHelper.gzFileToString(new File(out, "unmated.sam.gz"));
//      assertEquals(expected_unmated.replaceAll("^@.*$", ""), unmated.replaceAll("^@.*$", ""));
//
//      final String expected_mated = FileUtils.fileToString("Y:\\users\\mehul\\workspace\\rtg\\internal\\com\\rtg\\regression\\data\\long_pe_rammap\\bug997\\expect_out\\mwted.sam");
//      final String mated = FileHelper.gzFileToString(new File(out, "mated.sam.gz"));
//      assertEquals(expected_mated.replaceAll("^@.*$", ""), mated.replaceAll("^@.*$", ""));
//    } finally {
//      assertTrue(FileHelper.deleteAll(dir));
//    }
//  }

  public void testIonTorrentPairPens() throws Exception {
    try (final TestDirectory outer = new TestDirectory("map-end2end")) {
      final String templateStr = ">t" + LS + "GGTCCGTCACCTCCTGGTATTCCAGTTTGATGTTATCGGCAATACACAGCGCATTCATCAGTGGCTGGACGGTTTTTTCCTGCATATCCACGAGGTAGTAGTAACTGCCATCCTCGCAGCCATCGGGCGACTGGCGGGTACGTAAAACCTGCAAGGTGGGACCGGATAAGTAGTCAACCTGTGACCATTCTTCGCTGACATCATCCTGGTGGCTGACCACATCACTGAAACGCGGTGGGGTGAGGTCCTTAAATTTACTGATGGCCTTCAGGTCATCACTTCTGTCATCGCATGCGATGA" + LS;
      final String read1 = ">left" + LS + "GCCAGTCGCCCGATGGCTGCGAGGATGGCAGTTACTACTACCTCGTGGATATGCAGGAAAAAACCGTCCAGCCACTGATGAATGCGCTGTGTATTGCCGATAACATCAAGCTGAGGA" + LS;
      //    final String read1 = ">left" + StringUtils.LS + "TCCTCAGCTTGATGTTATCGGCAATACACAGCGCATTCATCAGTGGCTGGACGGTTTTTTCCTGCATATCCACGAGGTAGTAGTAACTGCCATCCTCGCAGCCATCGGGCGACTGGC" + StringUtils.LS;
      final String read2 = ">righ" + LS + "TTGATGTTATCGGCAATACACAGCGCATTCATCAGTGGCTGGACGGTTTTTTCCTGCATATCCACGAGGTAGTAGTAACTGCCATCCTCGCAGCCATCGGGCGACTGGCCTGAGTCGGAGACAC" + LS;

      final File template = new File(outer, "template");
      final File reads = new File(outer, "Reads");
      final File left = new File(reads, "left");
      final File right = new File(reads, "right");
      final File out = new File(outer, "out");

      ReaderTestUtils.getReaderDNA(templateStr, template, null).close();
      ReaderTestUtils.getReaderDNA(read1, left, new SdfId(123L)).close();
      ReaderTestUtils.getReaderDNA(read2, right, new SdfId(123L)).close();

      final MainResult r = MainResult.run(new MapCli(), "-i", reads.getPath(),
        "-t", template.getPath(),
        "-o", out.getPath(),
        "--" + MapFlags.MISMATCH_PENALTY_FLAG, "9",
        "--" + MapFlags.GAP_OPEN_PENALTY_FLAG, "19",
        "--" + MapFlags.GAP_EXTEND_PENALTY_FLAG, "1",
        "--" + MapFlags.ALIGNER_BAND_WIDTH_FACTOR_FLAG, "0.75",
        "--" + MapFlags.SOFT_CLIP_DISTANCE_FLAG, "0",
        "--" + MapFlags.SAM_FLAG,
        "--" + MapFlags.DONT_UNIFY_FLAG,
        "--" + MapFlags.ALIGNER_MODE_FLAG, "general");
      assertEquals(r.err(), 0, r.rc());

      final String mated = FileHelper.gzFileToString(new File(out, "mated.sam.gz"));
      TestUtils.containsAll(mated,
        "0\t83\tt\t26\t55\t1=7I109=\t=\t27\t110\tTCCTCAGCTTGATGTTATCGGCAATACACAGCGCATTCATCAGTGGCTGGACGGTTTTTTCCTGCATATCCACGAGGTAGTAGTAACTGCCATCCTCGCAGCCATCGGGCGACTGGC\t*\tAS:i:26\tNM:i:7\tXA:i:60\tIH:i:1\tNH:i:1\n",
        "0\t163\tt\t27\t55\t109=15I\t=\t26\t-110\tTTGATGTTATCGGCAATACACAGCGCATTCATCAGTGGCTGGACGGTTTTTTCCTGCATATCCACGAGGTAGTAGTAACTGCCATCCTCGCAGCCATCGGGCGACTGGCCTGAGTCGGAGACAC\t*\tAS:i:34\tNM:i:15\tXA:i:60\tIH:i:1\tNH:i:1\n"
      );
    }
  }

  public void testMappingAtSamePosition() throws Exception {
    try (TestDirectory td = new TestDirectory()) {
      final File template = new File(td, "tpl.fa.gz");
      FileHelper.resourceToGzFile("com/rtg/ngs/resources/ref_1_46877260.fasta", template);

      final File reads = new File(td, "reads.fa");
      FileHelper.resourceToFile("com/rtg/ngs/resources/1_46877260_xxx.sam", reads);

      final File fullSdf = new File(td, "sdf_full");
      final FormatCli format = new FormatCli();
      final MemoryPrintStream mps = new MemoryPrintStream();
      final int rc = format.mainInit(new String[]{"-o", fullSdf.getAbsolutePath(), template.getAbsolutePath()}, mps.outputStream(), mps.printStream());
      assertEquals(mps.toString(), 0, rc);

      final File out = new File(td, "out");

      final MapCli map = new MapCli();
      final String[] args = {
                                                "-i", reads.getPath(),
                                                "-t", fullSdf.getPath(),
                                                "-o", out.getPath(),
                                                "-F", "sam-pe",
                                                "--" + MapFlags.SAM_FLAG,
                                                "--" + MapFlags.DONT_UNIFY_FLAG
      };

      final int code = map.mainInit(args, mps.outputStream(), mps.printStream());
      assertEquals(mps.toString(), 0, code);

      final String mated = FileHelper.gzFileToString(new File(out, "mated.sam.gz"));
      final String actualFixed = TestUtils.stripSAMHeader(mated);

      mNano.check("matesamepos.sam", actualFixed, true);
    }
  }

  public void testTableLoading() throws Exception {
    try (final TestDirectory outer = new TestDirectory("map-end2end")) {
      final String templateStr = ">t" + LS + "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATATTTCTCCACATCCTCCCCAACACCTGTTGCTTCCTGACTTTTTAATGATCACCATTCTAACAGGTGTGAGATGGTATCTCACTGTGGTTTTGATTTGCATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" + LS;
      final String read1 = ">left" + LS + "TATTTCTCCACATCCTCCCCAACACCTGTTGCTTCCTGAATTTTAATGATTGCCATTCTAACTGGTGTGAGATGGTATCTCACTGTGGTTTTAATTTGCAT" + LS;

      final File template = new File(outer, "template");
      final File reads = new File(outer, "Reads");
      final File left = new File(reads, "left");
      final File out = new File(outer, "out");

      ReaderTestUtils.getReaderDNA(templateStr, template, null).close();
      ReaderTestUtils.getReaderDNA(read1, left, new SdfId(123L)).close();

      final File samrg = new File(outer, "header");
      FileUtils.stringToFile("@RG\tID:NA12878-1\tSM:NA12878\tPL:ILLUMINA", samrg);

      final MainResult r = MainResult.run(new MapCli(), "-i", left.getPath(),
        "-t", template.getPath(),
        "-o", out.getPath(),
        "--" + MapFlags.SAM_FLAG,
        "-Z",
        "--aligner-mode", "table",
        "--sam-rg", samrg.getPath());
      assertEquals(r.err(), 0, r.rc());

      final String mated = FileHelper.fileToString(new File(out, "alignments.sam"));

      //the default alignment penalties for non-illumina should make this alignment score 65 (20 for 1 long delete). The table specifies 1 long delete as 22, for a total of AS 67.
      //if any of these assumptions are false, this exact test probably won't make sense.
      assertEquals(20, EditDistanceFactory.DEFAULT_GAP_OPEN_PENALTY + EditDistanceFactory.DEFAULT_GAP_EXTEND_PENALTY);
      assertEquals("alignertable", EditDistanceFactory.DEFAULT_SINGLE_INDEL_TABLE);
      TestUtils.containsAll(FileHelper.resourceToString(PropertiesUtils.PropertyType.ALIGNMENT_PROPERTY_TYPE.path() + EditDistanceFactory.DEFAULT_SINGLE_INDEL_TABLE + ".properties"), "error_del_penalty=22,");

      TestUtils.containsAll(mated, "0\t0\tt\t32\t37\t39=1X4=1D6=2X10=1X29=1X8=\t*\t0\t0\tTATTTCTCCACATCCTCCCCAACACCTGTTGCTTCCTGAATTTTAATGATTGCCATTCTAACTGGTGTGAGATGGTATCTCACTGTGGTTTTAATTTGCAT\t*\tAS:i:67\tNM:i:6\tRG:Z:NA12878-1\tIH:i:1\tNH:i:1\n"
      );
    }
  }
}
