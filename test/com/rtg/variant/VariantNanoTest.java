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
import static com.rtg.sam.SharedSamConstants.SAM1_AMB;
import static com.rtg.sam.SharedSamConstants.SAM3;
import static com.rtg.sam.SharedSamConstants.SAM3_LENGTH;
import static com.rtg.sam.SharedSamConstants.SAM9;
import static com.rtg.sam.SharedSamConstants.SAM9_LENGTH;
import static com.rtg.sam.SharedSamConstants.SAMHEADER1;
import static com.rtg.sam.SharedSamConstants.SAM_LENGTH;
import static com.rtg.sam.SharedSamConstants.SAM_M;
import static com.rtg.sam.SharedSamConstants.SAM_M_LENGTH;
import static com.rtg.sam.SharedSamConstants.TEMPLATE_SDF_ID;
import static com.rtg.util.StringUtils.FS;
import static com.rtg.util.StringUtils.LS;
import static com.rtg.util.StringUtils.TAB;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractNanoTest;
import com.rtg.launcher.MainResult;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.sam.SamFilterOptions;
import com.rtg.sam.SharedSamConstants;
import com.rtg.tabix.TabixIndexer;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.Utils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.BgzipFileHelper;
import com.rtg.util.test.FileHelper;
import com.rtg.variant.bayes.multisample.AbstractMultisampleCli;
import com.rtg.variant.bayes.multisample.singleton.SingletonCli;

/**
 */
public class VariantNanoTest extends AbstractNanoTest {

  private static final String REF_SEQS67 = ""
      + ">g1" + LS
      +  SharedSamConstants.REF_BODY + LS
      + "agcatcgatcagcta\n"
      + ">gempty" + LS
      + LS
      ;

  public void test1() throws Exception {
    final String[] args0 = {"-a", "--Xno-complex-calls", "--snps-only"};
    check(REF_SEQS, SAM1, 1, null, 0, SAM_LENGTH, args0);
  }

  // output with default -T
  public void test2() throws Exception {
    final String[] args0 = {"--Xno-complex-calls", "--snps-only"};
    check(REF_SEQS, SAM1, 2, null, 0, SAM_LENGTH, args0);
  }

  public void test3() throws Exception {
    final String[] args0 = {"-a", "--Xno-complex-calls", "--snps-only"};
    check(REF_SEQS, SAM3, 3, null, 0, SAM3_LENGTH, args0);
  }

  // Include insertions and deletions
  // aa atcg actg gtca gcta gg
  // atcg actg
  // atcg g gtc deletion
  // cg actg tt
  // g actg ctc
  // g aTTT ctg insertion read
  // g actg insertion effect
  // g actg ctc
  // ttca gcta
  // e
  private static final String SAM4 = ""
      + "@HD" + TAB + "VN:1.0" + TAB + "SO:coordinate" + LS
      + TEMPLATE_SDF_ID + LS
      + "@SQ" + TAB + "SN:g1" + TAB + "LN:20" + LS
      + "@SQ" + TAB + "SN:gempty" + TAB + "LN:0" + LS
      + "0" + TAB + "0" + TAB + "g1" + TAB + "3" + TAB + "255" + TAB + "8M"     + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ATCGACTG" + TAB + "````````" + TAB + "AS:i:0" + LS
      + "1" + TAB + "0" + TAB + "g1" + TAB + "3" + TAB + "255" + TAB + "4M3D4M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ATCGGGTC" + TAB + "````````" + TAB + "AS:i:0" + LS
      + "2" + TAB + "0" + TAB + "g1" + TAB + "5" + TAB + "255" + TAB + "8M"     + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACTGTT" + TAB + "````````" + TAB + "AS:i:1" + LS
      + "3" + TAB + "0" + TAB + "g1" + TAB + "6" + TAB + "255" + TAB + "8M"     + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "GACTGCTC" + TAB + "````````" + TAB + "AS:i:1" + LS
      + "4" + TAB + "0" + TAB + "g1" + TAB + "6" + TAB + "255" + TAB + "2M3I3M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "GAtttCTG" + TAB + "````````" + TAB + "AS:i:1" + LS
      + "5" + TAB + "0" + TAB + "g1" + TAB + "11" + TAB + "255" + TAB + "8M"    + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "TTCAGCTA" + TAB + "````````" + TAB + "AS:i:1" + LS;

  public void test4() throws Exception {
    final String[] args0 = {"-a", "--Xno-complex-calls", "--snps-only"};
    check(REF_SEQS, SAM4, 4, null, 0, 6 * 8, args0);
  }

  // aa atcg actg gtca gcta gg
  // atcg actg
  // atcg actg
  // cg actg tt
  // g actg ctc
  // g actg ctc
  // ttca gcta
  // e
  // vary quality scores
  private static final String SAM5 = ""
      + "@HD" + TAB + "VN:1.0" + TAB + "SO:coordinate" + LS
      + TEMPLATE_SDF_ID + LS
      + "@SQ" + TAB + "SN:g1" + TAB + "LN:20" + LS
      + "@SQ" + TAB + "SN:gempty" + TAB + "LN:0" + LS
      + "0" + TAB + "0" + TAB + "g1" + TAB + "3"  + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ATCGACTG" + TAB + "&'(`````" + TAB + "AS:i:0" + LS
      + "1" + TAB + "0" + TAB + "g1" + TAB + "3"  + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ATCGACTG" + TAB + "````````" + TAB + "AS:i:0" + LS
      + "2" + TAB + "0" + TAB + "g1" + TAB + "5"  + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACTGTT" + TAB + "````````" + TAB + "AS:i:1" + LS
      + "3" + TAB + "0" + TAB + "g1" + TAB + "6"  + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "GACTGCTC" + TAB + "````````" + TAB + "AS:i:1" + LS
      + "4" + TAB + "0" + TAB + "g1" + TAB + "6"  + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "GACTGCTC" + TAB + "````````" + TAB + "AS:i:1" + LS
      + "5" + TAB + "0" + TAB + "g1" + TAB + "11" + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "TTCAGCTA" + TAB + "````````" + TAB + "AS:i:1" + LS;

  public void test5() throws Exception {
    final String[] args0 = {"-a", "--Xno-complex-calls", "--snps-only"};
    check(REF_SEQS, SAM5, 5, null, 0, 6 * 8, args0);
  }

  // aa atcg actg gtca gcta gg
  // atcg actg
  // atcg actg
  // cg actg tt
  // g actg ctc
  // g actg ctc
  // ttca gcta
  // e
  // no quality scores - use default
  private static final String SAM6 = ""
      + "@HD" + TAB + "VN:1.0" + TAB + "SO:coordinate" + LS
      + TEMPLATE_SDF_ID + LS
      + "@SQ" + TAB + "SN:g1" + TAB + "LN:20" + LS
      + "@SQ" + TAB + "SN:gempty" + TAB + "LN:0" + LS
      + "0" + TAB + "0" + TAB + "g1" + TAB + "3"  + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ATCGACTG" + TAB + "*" + TAB + "AS:i:0" + LS
      + "1" + TAB + "0" + TAB + "g1" + TAB + "3"  + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ATCGACTG" + TAB + "*" + TAB + "AS:i:0" + LS
      + "2" + TAB + "0" + TAB + "g1" + TAB + "5"  + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACTGTT" + TAB + "*" + TAB + "AS:i:1" + LS
      + "3" + TAB + "0" + TAB + "g1" + TAB + "6"  + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "GACTGCTC" + TAB + "*" + TAB + "AS:i:1" + LS
      + "4" + TAB + "0" + TAB + "g1" + TAB + "6"  + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "GACTGCTC" + TAB + "*" + TAB + "AS:i:1" + LS
      + "5" + TAB + "0" + TAB + "g1" + TAB + "11" + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "TTCAGCTA" + TAB + "*" + TAB + "AS:i:1" + LS;
  private static final int SAM6_LENGTH = 6 * 8;

  public void test6() throws Exception {
    final String[] args0 = {"-a", "--Xno-complex-calls", "--snps-only"};
    check(REF_SEQS, SAM6, 6, null, 0, SAM6_LENGTH, args0);
  }

  public void test7() throws Exception {
    final String[] args0 = {"-a", "-q", "10", "--Xno-complex-calls", "--snps-only"};
    check(REF_SEQS, SAM6, 7, null, 0, SAM6_LENGTH, args0);
  }

  // multiple sequences
  public void test8() throws Exception {
    final String[] args0 = {"-a", "--Xno-complex-calls", "--snps-only"};
    check(REF_SEQS_M, SAM_M, 8, null, 0, SAM_M_LENGTH, args0);
  }

  //  // multiple sequences will move to <code>AbstractMultisampleCli</code>
  //  public void testUnsorted() {
  //    final String[] args0 = {"-a"};
  //    check(REF_SEQS, SAM_UNSORTED, args0, 8, "sort", "", 1, 10, false, null, false);
  //  }

  public void test9() throws Exception {
    final String[] args0 = {"-a", "--filter-ambiguity", "100%", "--Xno-complex-calls", "--snps-only"};
    check(REF_SEQS, SAM9, 9, null, 0, SAM9_LENGTH, args0);
  }

  public void testIHFilter() throws Exception {
    final String[] args0 = {"-a", "--" + SamFilterOptions.MAX_HITS_FLAG, "2", };
    check(REF_SEQS, SAM9, -1, null, "2 records skipped due to input filtering criteria", 0, 10, null, true, SAM9_LENGTH, args0);
  }

  public void testASFilter() throws Exception {
    final String[] args0 = {"-a", "--max-as-unmated", "0"};
    check(REF_SEQS, SAM9, -1, null, "4 records skipped due to input filtering criteria", 0, 10, null, true, SAM9_LENGTH, args0);
  }

  private void check(final String refSeq, final String sam, final int expNum, final String errorMsg, final int errCode, final long usageExp, final String... args0) throws Exception {
    check(refSeq, sam, expNum, errorMsg, "", errCode, 10, null, true, usageExp, args0);
  }

  private void check(String refSeq, String sam, int expNum, String errorMsg, int errCode, int intSeparation, String refFile, long usageExp, String... args0) throws Exception {
    check(refSeq, sam, expNum, errorMsg, "", errCode, intSeparation, refFile, true, usageExp, args0);
  }

  private void check(String refSeq, String sam, int expNum, String errorMsg, String logMsg, int errCode, int intSeparation, String refFile, boolean index, long usageExp, String... args0) throws Exception {
    try (final TestDirectory testdir = new TestDirectory("variance_check")) {
      final File output = new File(testdir, "variant_out");
      final File input = new File(testdir, "variant_in");
      if (!input.mkdir()) {
        throw new IOException("Whut, couldn't make test directory");
      }
      final File templ = new File(testdir, "variant_template");
      if (!templ.mkdir()) {
        throw new IOException("Whut, couldn't make test directory");
      }
      // FileUtils.saveFile(new File(input, OUT_SAM), sam);
      final File file = new File(input, OUT_SAM + ".gz");
      BgzipFileHelper.streamToBgzipFile(new ByteArrayInputStream(sam.getBytes()), file);
      if (index) {
        new TabixIndexer(file, new File(input, OUT_SAM + ".gz.tbi")).saveSamIndex();
      }
      final String outn = output.getPath();
      final String inn = input.getPath();
      ReaderTestUtils.getDNADir(refSeq, templ);
      if (refFile != null) {
        FileUtils.stringToFile(refFile, new File(templ, "reference.txt"));
      }
      final String[] alignArgs = {"-t", templ.getPath(),
        "-o", outn, "--Xpriors", "testhumanprior", "-Z", "--Xinteresting-separation",
        Integer.toString(intSeparation), inn + FS + OUT_SAM + ".gz", "--" + AbstractMultisampleCli.NO_CALIBRATION
      };
      final String[] args1;
      if (args0 != null) {
        args1 = Utils.append(alignArgs, args0);
      } else {
        args1 = alignArgs;
      }
      final String[] argsv;
      argsv = args1;
      final String[] args;
      if (!Arrays.asList(argsv).contains("-m")) {
        args = Utils.append(argsv, "-m", "default");
      } else {
        args = argsv;
      }
      final MainResult r = MainResult.run(new SingletonCli(), args);
      assertEquals(r.err(), errCode, r.rc());
      // assertEquals("", out.toString()); Now outputting stats
      //          System.err.println(errStr);
      if (errorMsg != null) {
        TestUtils.containsAll(r.err(), errorMsg);
      } else {
        assertEquals(r.err(), 0, r.err().length());
      }
      if (logMsg != null && !logMsg.equals("")) {
        final String lout = FileUtils.fileToString(new File(outn, "snp.log"));
        assertTrue(lout, lout.contains(logMsg));
      }
      final String usageLog = AbstractCli.lastUsageLog();
      //System.err.println(usageLog);
      if (usageExp >= 0) {
        TestUtils.containsAll(usageLog, "[Usage beginning module=snp runId=", ", Usage end module=snp runId=", " metric=" + usageExp + " success=true]");
      }
      if ((errorMsg == null || "".equals(errorMsg)) && (expNum >= 0)) {
        final String result = FileUtils.fileToString(new File(output, VariantParams.VCF_OUT_SUFFIX));
        final String actualFixed = TestUtils.sanitizeVcfHeader(result);
        //System.err.println(actual);
        //final String expectedFinal = FileHelper.resourceToString("com/rtg/variant/resources/variantnanotest" + expNum + ".vcf").replaceAll("\r", "");
        //final String actualFinal = actualFixed.replaceAll("\r", "");
        //System.out.println("Expected: \n" + expectedFinal);
        //System.out.println("Actual: \n" + actual);
        //assertEquals(expectedFinal, actualFinal);
        mNano.check("variantnanotest" + expNum + ".vcf", actualFixed, false);
      }
    }
  }

  private void checkMultifile(final String refSeq, final String[] args0, final int expNum, final String maps1, final String maps2, final String maps3) throws Exception {
    try (final TestDirectory testdir = new TestDirectory("variance_check")) {
      final File output = new File(testdir, "variant_out");
      final File input = new File(testdir, "variant_in");
      if (!input.mkdir()) {
        throw new IOException("Whut, couldn't make test directory");
      }
      final File templ = new File(testdir, "variant_template");
      if (!templ.mkdir()) {
        throw new IOException("Whut, couldn't make test directory");
      }
      final File sam1 = new File(input, "sam1.sam.gz");
      VariantTestUtils.bgzipAndIndex(maps1, sam1);
      final File sam2 = new File(input, "sam2.sam.gz");
      VariantTestUtils.bgzipAndIndex(maps2, sam2);
      final File sam3 = new File(input, "sam3.sam.gz");
      VariantTestUtils.bgzipAndIndex(maps3, sam3);

      final String outn = output.getPath();
      ReaderTestUtils.getDNADir(refSeq, templ);
      final String[] alignArgs = {"-t", templ.getPath(), "-o", outn, "--Xpriors", "testhumanprior", "-m", "default", "--keep-duplicates", "-Z", sam1.getPath(), sam2.getPath(), sam3.getPath(), "--" + AbstractMultisampleCli.NO_CALIBRATION};
      final String[] args = Utils.append(alignArgs, args0);
      final MainResult r = MainResult.run(new SingletonCli(), args);
      assertEquals(r.err(), 0, r.rc());

      final String result = FileUtils.fileToString(new File(output, VariantParams.VCF_OUT_SUFFIX));
      final String actualFixed = TestUtils.sanitizeVcfHeader(result);
      mNano.check("variantnanotest" + expNum + ".vcf", actualFixed, false);
    }
  }

  public void test5Compressed() throws Exception {
    final String[] args0 = {"-a", "--Xno-complex-calls", "--snps-only"};
    checkCompressed(SAM5, args0, 5);
  }

  private void checkCompressed(final String sam, final String[] args0, final int expNum) throws Exception {
    try (final TestDirectory dir = new TestDirectory("testcheckCompressed")) {
      final File output = new File(dir, "variant_out");
      final File file = new File(dir, OUT_SAM + ".gz");
      BgzipFileHelper.streamToBgzipFile(new ByteArrayInputStream(sam.getBytes()), file);
      new TabixIndexer(file, new File(dir, OUT_SAM + ".gz.tbi")).saveSamIndex();
      final String outn = output.getPath();
      final String inn = dir.getPath();
      final File templ = ReaderTestUtils.getDNADir(REF_SEQS);
      try {
        final String[] alignArgs = {"-t", templ.getPath(), "-o", outn, "--Xpriors", "testhumanprior", "-m", "default", "-Z", "--keep-duplicates", inn + FS + OUT_SAM + FileUtils.GZ_SUFFIX, "--" + AbstractMultisampleCli.NO_CALIBRATION};
        final String[] args = Utils.append(alignArgs, args0);
        final MainResult r = MainResult.run(new SingletonCli(), args);
        assertEquals("", r.err());
        final String result = FileUtils.fileToString(new File(output, VariantParams.VCF_OUT_SUFFIX));
        final String actual = result.replace("Version", "");
        //System.err.println(actual);
        assertTrue(actual.startsWith("##fileformat="));
        final String actualFixed = TestUtils.sanitizeVcfHeader(actual);
        mNano.check("variantnanotest" + expNum + ".vcf", actualFixed, false);
      } finally {
        FileHelper.deleteAll(templ);
      }
    }
  }

  static final String SAM_HEAD1 = "" + "@HD" + TAB + "VN:1.0" + TAB + "SO:coordinate\n" + "@SQ" + TAB + "SN:gi" + TAB + "LN:30\n";

  static final String SAM_REC_OK1 = "" + "okread" + TAB + "0" + TAB + "gi" + TAB + "2" + TAB + "255" + TAB + "10M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "AAAAAAAAAA" + TAB + "IB7?*III<I" + TAB + "AS:i:0" + TAB + "IH:i:1" + LS;
  // Invalid IH value
  static final String SAM_REC_BAD1 = "" + "okread" + TAB + "0" + TAB + "gi" + TAB + "2" + TAB + "255" + TAB + "10M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "AAAAAAAAAA" + TAB + "IB7?*III<I" + TAB + "AS:i:0" + TAB + "IH:i:0" + LS;
  // bad syntax in field
  static final String SAM_REC_BAD2 = "" + "67" + TAB + "16" + TAB + "gi" + TAB + "3" + TAB + "255" + TAB + "10M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "AAAAAAAAAA" + TAB + "GC@=I3IIII" + TAB + "AS:i:0" + TAB + "IH:0" + LS;
  // missing field
  static final String SAM_REC_BAD_FATAL = "" + "badread3a" + TAB + "16" + TAB + "gi0";

  static final String SAM_REC_BAD_CIGAR = "" + "badcigar" + TAB + "16" + TAB + "gi" + TAB + "3" + TAB + "255" + TAB + "2M1D1D8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "AAAAAAAAAA" + TAB + "GC@=I3IIII" + TAB + "AS:i:0" + TAB + "IH:i:1" + LS;

  static final String SAM_REC_GOOD_CIGAR = "" + "goodcigar" + TAB + "16" + TAB + "gi" + TAB + "3" + TAB + "255" + TAB + "2M1D1N7M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "AAAAAAAAAA" + TAB + "GC@=I3IIII" + TAB + "AS:i:0" + TAB + "IH:i:1" + LS;


  static final String SAM_REC_SOFTCLIPPED = "" + "23" + TAB + "0" + TAB + "gi" + TAB + "23" + TAB + "255" + TAB + "8M2S" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "AAAAAAAAAA" + TAB + "IB7?*III<I" + TAB + "AS:i:0" + TAB + "IH:i:1" + LS;

  static final String SAM_REC_SOFTCLIPPED2 = "" + "23" + TAB + "0" + TAB + "gi" + TAB + "1" + TAB + "255" + TAB + "2S8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "AAAAAAAAAA" + TAB + "IB7?*III<I" + TAB + "AS:i:0" + TAB + "IH:i:1" + LS;

  static final String SAM_REC_NEGATIVEPOS = "" + "23" + TAB + "0" + TAB + "gi" + TAB + "-1" + TAB + "255" + TAB + "2S8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "AAAAAAAAAA" + TAB + "IB7?*III<I" + TAB + "AS:i:0" + TAB + "IH:i:1" + LS;
  static final String SAM_REC_ZERO_POS = "" + "23" + TAB + "0" + TAB + "gi" + TAB + "0" + TAB + "255" + TAB + "2S8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "AAAAAAAAAA" + TAB + "IB7?*III<I" + TAB + "AS:i:0" + TAB + "IH:i:1" + LS;

  private static final String REF_DOESNTMATTER = "" + ">gi" + LS
      // 1234567890
      // 123456789012345678901234567890
      //
      + "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaa" + LS;
  // 123456789012345678901234567890

  private static final String REF_WRONGNAME = "" + ">g111" + LS + "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaa" + LS;

  private void checkSamFormatNoError(final String refSeq, final String sam, final String[] warnings, final String[] logs) throws Exception {
    try (final TestDirectory dir = new TestDirectory("testchecksamformat")) {
      final File output = new File(dir, "variant_out");
      // FileUtils.saveFile(new File(alignments, OUT_SAM), sam);
      final File file = new File(dir, OUT_SAM + ".gz");
      BgzipFileHelper.streamToBgzipFile(new ByteArrayInputStream(sam.getBytes()), file);
      new TabixIndexer(file, new File(dir, OUT_SAM + ".gz.tbi")).saveSamIndex();
      final String align = dir.getPath();
      final String outn = output.getPath();
      final File templ = ReaderTestUtils.getDNADir(refSeq);
      try {
        final MainResult r = MainResult.run(new SingletonCli(), "-t", templ.getPath(), "-m", "default", "-o", outn, align + FS + OUT_SAM + ".gz", "--" + AbstractMultisampleCli.NO_CALIBRATION);
        final String log = FileUtils.fileToString(new File(output, "snp.log"));
        TestUtils.containsAll(log, logs);
        if (warnings != null) {
          TestUtils.containsAll(r.err(), warnings);
        }
      } finally {
        assertTrue(FileHelper.deleteAll(templ));
      }
    }
  }

  private void checkSamFormatFatal(final String refSeq, final String sam, final String... errExpected) throws Exception {
    try (final TestDirectory dir = new TestDirectory("testchecksamformat")) {
      final File output = new File(dir, "variant_out");
      FileUtils.stringToFile(sam, new File(dir, OUT_SAM));
      final String align = dir.getPath();
      final String outn = output.getPath();
      final File templ = ReaderTestUtils.getDNADir(refSeq);
      try {
        final MainResult r = MainResult.run(new SingletonCli(), "-t", templ.getPath(), "-o", outn, "-m", "illumina", align + FS + OUT_SAM, "--" + AbstractMultisampleCli.NO_CALIBRATION);
        assertEquals(1, r.rc());
        TestUtils.containsAll(r.err(), errExpected);
        assertFalse(r.err().contains("RTG has encountered a diff"));
      } finally {
        FileHelper.deleteAll(templ);
      }
    }
  }

  public void testSamFormatOk() throws Exception {
    // no warning
    final String sam = SAM_HEAD1 + SAM_REC_OK1;
    checkSamFormatNoError(REF_DOESNTMATTER, sam, null, new String[] {});
  }

  public void testSamFormatOk2() throws Exception {
    // no warning
    final String sam = SAM_HEAD1 + SAM_REC_OK1 + SAM_REC_OK1 + SAM_REC_GOOD_CIGAR;
    checkSamFormatNoError(REF_DOESNTMATTER, sam, null, new String[] {});
  }

  private String[] makeSkippedWarning(final int num, final int allNum) {
    final String skipped = " records skipped because of SAM format problems.";
    final List<String> list = new ArrayList<>();
    if (num > 0) {
      list.add(num + skipped);
    }
    if (allNum > 0) {
      list.add(allNum + skipped);
    }
    return list.toArray(new String[list.size()]);
  }

  public void testSamFormatBadIH() throws Exception {
    final String sam = SAM_HEAD1 + SAM_REC_OK1 + SAM_REC_BAD1;
    checkSamFormatNoError(REF_DOESNTMATTER, sam, makeSkippedWarning(0, 0), new String[] {"1 records skipped due to input filtering criteria"});
  }

  public void testSamFormatFirstRecBad() throws Exception {
    final String sam = SAM_HEAD1 + SAM_REC_BAD2;
    checkSamFormatNoError(REF_DOESNTMATTER, sam, makeSkippedWarning(1, 0), new String[] {});
  }

  public void testSamFormatSecondRecBad() throws Exception {
    final String sam = SAM_HEAD1 + SAM_REC_OK1 + SAM_REC_BAD2;
    checkSamFormatNoError(REF_DOESNTMATTER, sam, makeSkippedWarning(1, 1), new String[] {});
  }

  public void testSamFormatThirdRecBad() throws Exception {
    final String sam = SAM_HEAD1 + SAM_REC_OK1 + SAM_REC_OK1 + SAM_REC_BAD2;
    checkSamFormatNoError(REF_DOESNTMATTER, sam, makeSkippedWarning(1, 1), new String[] {});
  }

  public void testSamFormatOneRecBadFatal() throws Exception {
    final String sam = SAM_HEAD1 + SAM_REC_OK1 + SAM_REC_BAD_FATAL;
    checkSamFormatFatal(REF_DOESNTMATTER, sam, "SAM record has an irrecoverable problem in file ");
  }

  // static final String ERROR_CIGAR_IDN = ""
  // + "No M operator between pair of IDN operators in CIGAR";
  // no error at all
  public void testSamFormatBadCigarFirst() throws Exception {
    final String sam = SAM_HEAD1 + SAM_REC_BAD_CIGAR;
    checkSamFormatNoError(REF_DOESNTMATTER, sam, null, new String[] {}); // makeSkippedWarning(1,
    // 0));
    // "records skipped in file \"All Files\" because of SAM format problems");
  }

  public void testSamFormatBadCigarSecond() throws Exception {
    final String sam = SAM_HEAD1 + SAM_REC_OK1 + SAM_REC_BAD_CIGAR;
    checkSamFormatNoError(REF_DOESNTMATTER, sam, null, new String[] {}); // makeSkippedWarning(1,
    // 1));
  }

  public void testSamFormatOneBadCigarThird() throws Exception {
    final String sam = SAM_HEAD1 + SAM_REC_OK1 + SAM_REC_OK1 + SAM_REC_BAD_CIGAR;
    checkSamFormatNoError(REF_DOESNTMATTER, sam, null, new String[] {});
  }

  public void testSamFormatWrongSequ() throws Exception {
    final String sam = SAM_HEAD1 + SAM_REC_OK1;
    checkSamFormatFatal(REF_WRONGNAME, sam, "Sequence 'gi' not in reference SDF");
  }

  public void testSoftClippedRight() throws Exception {
    final String sam = SAM_HEAD1 + SAM_REC_SOFTCLIPPED;
    checkSamFormatNoError(REF_DOESNTMATTER, sam, null, new String[] {});
  }

  public void testSoftClippedLeft() throws Exception {
    final String sam = SAM_HEAD1 + SAM_REC_SOFTCLIPPED2;
    checkSamFormatNoError(REF_DOESNTMATTER, sam, null, new String[] {});
  }

  public void testNegativePos() throws Exception {
    final String sam = SAM_HEAD1 + SAM_REC_NEGATIVEPOS;
    checkSamFormatNoError(REF_DOESNTMATTER, sam, makeSkippedWarning(0, 0), new String[] {"1 records skipped due to input filtering criteria"});
  }

  public void testZeroPos() throws Exception {
    final String sam = SAM_HEAD1 + SAM_REC_ZERO_POS;
    checkSamFormatNoError(REF_DOESNTMATTER, sam, makeSkippedWarning(0, 1), new String[] {});
  }

  // run same tests with records split up into 3 files
  public void test3Multifile() throws Exception {
    final String[] args0 = {"-a", "--Xno-complex-calls", "--snps-only"};
    checkMultifile(REF_SEQS, args0, 3, IN3SAM3HEADER + OK0 + OK1, IN3SAM3HEADER + OK2 + OK3, IN3SAM3HEADER + OK4 + OK5);
    checkMultifile(REF_SEQS, args0, 3, IN3SAM3HEADER + OK0 + OK4, IN3SAM3HEADER + OK1 + OK3, IN3SAM3HEADER + OK2 + OK5);
  }

  // check case when there are Ns in the genome
  private static final String REF_SEQS_N = "" + ">g1" + LS + "aa" + "atcg" + "actg" + "gNca" + "gcta" + "gg" + LS + ">gempty" + LS + LS;

  public void test10() throws Exception {
    final String[] args0 = {"-a", "--Xno-complex-calls", "--snps-only"};
    check(REF_SEQS_N, SAM1, 10, null, 0, SAM_LENGTH, args0);
  }

  public void test11() throws Exception {
    final String[] args0 = {"-a", "--Xno-complex-calls", "--snps-only"};
    check(SharedSamConstants.REF_SEQS11, SharedSamConstants.SAM11, 11, null, 0, SharedSamConstants.SAM_LENGTH11, args0);
  }

  // 12 3456 7890 1234 5678 90
  // aa atcg actg gtca gcta gg
  // atcg acdg
  // atcg acdg
  // cg acdg tt
  // g acdg ctc
  // g acdg ctc
  // ttca gcta
  // SAM format file after mapping and alignment
  private static final String SAM_BODY_D = ""
      + "0" + TAB + "0" + TAB + "g1" + TAB + "3"  + TAB + "255" + TAB + "6M1D1M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ATCGACG" + TAB + "```````" + TAB + "AS:i:0" + LS
      + "1" + TAB + "0" + TAB + "g1" + TAB + "3"  + TAB + "255" + TAB + "6M1D1M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ATCGACG" + TAB + "```````" + TAB + "AS:i:0" + LS
      + "2" + TAB + "0" + TAB + "g1" + TAB + "5"  + TAB + "255" + TAB + "4M1D3M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACGTT" + TAB + "```````" + TAB + "AS:i:1" + LS
      + "3" + TAB + "0" + TAB + "g1" + TAB + "6"  + TAB + "255" + TAB + "3M1D4M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "GACGCTC" + TAB + "```````" + TAB + "AS:i:1" + LS
      + "4" + TAB + "0" + TAB + "g1" + TAB + "6"  + TAB + "255" + TAB + "3M1D4M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "GACGCTC" + TAB + "```````" + TAB + "AS:i:1" + LS
      + "5" + TAB + "0" + TAB + "g1" + TAB + "11" + TAB + "255" + TAB + "8M"     + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "TTCAGCTA" + TAB + "````````" + TAB + "AS:i:1" + LS;
  private static final String SAM1_D = ""
      + "@HD" + TAB + "VN:1.0" + TAB + "SO:coordinate" + LS
      + "@SQ" + TAB + "SN:g1" + TAB + "LN:20" + LS
      + "@SQ" + TAB + "SN:gempty" + TAB + "LN:0" + LS
      + TEMPLATE_SDF_ID + LS
      + SAM_BODY_D;
  private static final int SAM1_D_LENGTH = 5 * 7 + 8;

  // homozygous delete
  public void test12() throws Exception {
    final String[] args0 = {"-a"};
    check(REF_SEQS, SAM1_D, 12, null, 0, SAM1_D_LENGTH, args0);
  }

  // 12345678 901234567890
  // aaatcgac**tggtcagctagg
  // atcgacgctg
  // atcgacgctg
  // cgacgctggt
  // gacgctggtc
  // gacgctggtc
  // ttcagcta
  // SAM format file after mapping and alignment
  private static final String SAM_BODY_I = ""
      + "0" + TAB + "0" + TAB + "g1" + TAB + "3"  + TAB + "255" + TAB + "6M2I2M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ATCGACGCTG" + TAB + "``````````" + TAB + "AS:i:0" + LS
      + "1" + TAB + "0" + TAB + "g1" + TAB + "3"  + TAB + "255" + TAB + "6M2I2M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ATCGACGCTG" + TAB + "``````````" + TAB + "AS:i:0" + LS
      + "2" + TAB + "0" + TAB + "g1" + TAB + "5"  + TAB + "255" + TAB + "4M2I4M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACGCTGGT" + TAB + "``````````" + TAB + "AS:i:1" + LS
      + "3" + TAB + "0" + TAB + "g1" + TAB + "6"  + TAB + "255" + TAB + "3M2I5M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "GACGCTGGTC" + TAB + "``````````" + TAB + "AS:i:1" + LS
      + "4" + TAB + "0" + TAB + "g1" + TAB + "6"  + TAB + "255" + TAB + "3M2I5M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "GACGCTGGTC" + TAB + "``````````" + TAB + "AS:i:1" + LS
      + "5" + TAB + "0" + TAB + "g1" + TAB + "11" + TAB + "255" + TAB + "8M"     + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "TTCAGCTA"   + TAB + "````````" + TAB + "AS:i:1" + LS;
  private static final long SAM1_I_LENGTH = 5 * 10 + 8;

  private static final String SAM1_I = ""
      + "@HD" + TAB + "VN:1.0" + TAB + "SO:coordinate" + LS
      + "@SQ" + TAB + "SN:g1" + TAB + "LN:20" + LS
      + "@SQ" + TAB + "SN:gempty" + TAB + "LN:0" + LS
      + TEMPLATE_SDF_ID + LS
      + SAM_BODY_I;

  // homozygous delete - now correct yay!
  public void test13() throws Exception {
    final String[] args0 = {"-a"};
    check(REF_SEQS, SAM1_I, 13, null, 0, SAM1_I_LENGTH, args0);
  }

  // homozygous delete no -a
  public void test14() throws Exception {
    check(REF_SEQS, SAM1_I, 14, null, 0, SAM1_I_LENGTH);
  }

  public void test15() throws Exception {
    check(REF_SEQS, SAM1_I, 15, null, 0, SAM1_I_LENGTH, "-a");
  }

  public void test16() throws Exception {
    check(REF_SEQS, SAM1_I, 16, null, 0, 5, null, SAM1_I_LENGTH, "-a");
  }

  // 12345678  90 1234567890
  // aaatcgac**tg*gtcagctagg
  // atcgacgctg
  // atcgacgctg
  // cgacgctg*gt
  // gacgctg*gtc
  // gacgctg*gtc
  // tgagtcag
  // tg*ttcagct
  // agtcagct
  // agtcagct
  // agtcagct
  // agtcagct
  // agtcagct
  // ttcagcta
  // ttcagcta
  // ttcagcta
  // ttcagcta
  // SAM format file after mapping and alignment
  private static final String SAM_BODY_I2 = ""
      + "0"  + TAB + "0" + TAB + "g1" + TAB + "3"  + TAB + "255" + TAB + "6M2I2M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ATCGACGCTG" + TAB + "``````````" + TAB + "AS:i:0" + LS
      + "1"  + TAB + "0" + TAB + "g1" + TAB + "3"  + TAB + "255" + TAB + "6M2I2M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ATCGACGCTG" + TAB + "``````````" + TAB + "AS:i:0" + LS
      + "2"  + TAB + "0" + TAB + "g1" + TAB + "5"  + TAB + "255" + TAB + "4M2I4M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACGCTGGT" + TAB + "``````````" + TAB + "AS:i:1" + LS
      + "3"  + TAB + "0" + TAB + "g1" + TAB + "6"  + TAB + "255" + TAB + "3M2I5M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "GACGCTGGTC" + TAB + "``````````" + TAB + "AS:i:1" + LS
      + "4"  + TAB + "0" + TAB + "g1" + TAB + "6"  + TAB + "255" + TAB + "3M2I5M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "GACGCTGGTC" + TAB + "``````````" + TAB + "AS:i:1" + LS
      + "6"  + TAB + "0" + TAB + "g1" + TAB + "9"  + TAB + "255" + TAB + "2M1I5M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "TGAGTCAG"   + TAB + "````````"   + TAB + "AS:i:2" + LS
      + "7"  + TAB + "0" + TAB + "g1" + TAB + "9"  + TAB + "255" + TAB + "8M"     + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "TGTTCAGC"   + TAB + "````````"   + TAB + "AS:i:2" + LS
      + "5"  + TAB + "0" + TAB + "g1" + TAB + "11" + TAB + "255" + TAB + "8M"     + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "TTCAGCTA"   + TAB + "````````"   + TAB + "AS:i:1" + LS
      + "8"  + TAB + "0" + TAB + "g1" + TAB + "11" + TAB + "255" + TAB + "1I7M"   + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "AGTCAGCT"   + TAB + "````````"   + TAB + "AS:i:2" + LS
      + "9"  + TAB + "0" + TAB + "g1" + TAB + "11" + TAB + "255" + TAB + "1I7M"   + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "AGTCAGCT"   + TAB + "````````"   + TAB + "AS:i:3" + LS
      + "10" + TAB + "0" + TAB + "g1" + TAB + "11" + TAB + "255" + TAB + "1I7M"   + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "AGTCAGCT"   + TAB + "````````"   + TAB + "AS:i:3" + LS
      + "11" + TAB + "0" + TAB + "g1" + TAB + "11" + TAB + "255" + TAB + "1I7M"   + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "AGTCAGCT"   + TAB + "````````"   + TAB + "AS:i:3" + LS
      + "12" + TAB + "0" + TAB + "g1" + TAB + "11" + TAB + "255" + TAB + "1I7M"   + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "AGTCAGCT"   + TAB + "````````"   + TAB + "AS:i:3" + LS
      + "13" + TAB + "0" + TAB + "g1" + TAB + "11" + TAB + "255" + TAB + "8M"     + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "TTCAGCTA"   + TAB + "````````"   + TAB + "AS:i:1" + LS
      + "14" + TAB + "0" + TAB + "g1" + TAB + "11" + TAB + "255" + TAB + "8M"     + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "TTCAGCTA"   + TAB + "````````"   + TAB + "AS:i:1" + LS
      + "15" + TAB + "0" + TAB + "g1" + TAB + "11" + TAB + "255" + TAB + "8M"     + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "TTCAGCTA"   + TAB + "````````"   + TAB + "AS:i:1" + LS;
  private static final long SAM1_I2_LENGTH = 5 * 10 + 11 * 8;

  private static final String SAM1_I2 = ""
      + "@HD" + TAB + "VN:1.0" + TAB + "SO:coordinate" + LS
      + TEMPLATE_SDF_ID + LS
      + "@SQ" + TAB + "SN:g1" + TAB + "LN:20" + LS
      + "@SQ" + TAB + "SN:gempty" + TAB + "LN:0" + LS
      + SAM_BODY_I2;

  public void test17default() throws Exception {
    final String[] args0 = {"-a", "--sex=either"};
    check(REF_SEQS, SAM1_I2, 17, null, "Processing g1", 0, 1/*
     * interesting
     * separation
     */, null, true, SAM1_I2_LENGTH, args0);
  }

  public void test17diploid() throws Exception {
    check(REF_SEQS, SAM1_I2, 17, "", "Processing g1", 0, 1, SharedSamConstants.REF_DIPLOID, true, SAM1_I2_LENGTH, "-a", "--sex=either");
  }

  private static final String SAM1_II = ""
      + "@HD" + TAB + "VN:1.0" + TAB + "SO:coordinate" + LS
      + "@SQ" + TAB + "SN:g1" + TAB + "LN:35" + LS
      + "@SQ" + TAB + "SN:gempty" + TAB + "LN:0" + LS
      + TEMPLATE_SDF_ID + LS
      + SAM_BODY_I;
  // TODO test explicitly setting sex and haploid/diploid on different
  // chromosomes

  public void test18() throws Exception {
    final String[] args0 = {"-a"};
    check(REF_SEQS67, SAM1_II, 18, null, 0, SAM1_I_LENGTH, args0);
  }

  // homozygous delete no -a
  // see test 14
  public void test19() throws Exception {
    final String[] args0 = {};
    check(REF_SEQS, SAM1_I, 19, null, 0, SAM1_I_LENGTH, args0);
  }

  // homozygous delete but don't attempt a complex call - see test13
  public void test20() throws Exception {
    final String[] args0 = {"-a", "--Xno-complex-calls"};
    // remove g1  9 x i x 0.0 5 0.050 GC  5 0.050 - in case in singleton caller
    check(REF_SEQS, SAM1_I, 20, null, 0, SAM1_I_LENGTH, args0);
  }

  // see test12 - add nonidentity posterior
  public void test21() throws Exception {
    final String[] args0 = {"-a"};
    check(REF_SEQS, SAM1_D, 21, null, 0, SAM1_D_LENGTH, args0);
  }

  // Use alternate priors - see test 1
  public void test23() throws Exception {
    final String[] args0 = {"-a", "--Xno-complex-calls", "--snps-only"};
    check(REF_SEQS, SAM1, 23, null, 0, SAM_LENGTH, args0);
  }

  /** This tests that we can set read mapping qualities via the command line */
  public void test25() throws Exception {
    final String[] args0 = {"-a",
        "--rdefault-mated", "10",
        "--rdefault-unmated", "10" // -R is enough because SAM6 is unmated
        , "--Xno-complex-calls", "--snps-only"};
    check(REF_SEQS, SAM6, 25, null, 0, SAM6_LENGTH, args0);
  }

  // 12 3456 7890 1234 5678 90
  // aa atcg actg gtca gcta gg
  // atcN Nctg
  // atcN Nctg
  // cg aNNg tt
  // g acNN ctc
  // g acNN ctc
  // ttcN Ncta

  /** SAM format file after mapping and alignment **/
  public static final String SAM_BODY_26 = ""
      + "0" + TAB + "0" + TAB + "g1" + TAB + "3"  + TAB + "255" + TAB + "3M2N3M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ATCCTG" + TAB + "``````" + TAB + "AS:i:0" + LS
      + "1" + TAB + "0" + TAB + "g1" + TAB + "3"  + TAB + "255" + TAB + "3M2N3M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ATCCTG" + TAB + "``````" + TAB + "AS:i:0" + LS
      + "2" + TAB + "0" + TAB + "g1" + TAB + "5"  + TAB + "255" + TAB + "3M2N3M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGAGTT" + TAB + "``````" + TAB + "AS:i:1" + LS
      + "3" + TAB + "0" + TAB + "g1" + TAB + "6"  + TAB + "255" + TAB + "3M2N3M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "GACCTC" + TAB + "``````" + TAB + "AS:i:1" + LS
      + "4" + TAB + "0" + TAB + "g1" + TAB + "6"  + TAB + "255" + TAB + "3M2N3M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "GACCTC" + TAB + "``````" + TAB + "AS:i:1" + LS
      + "5" + TAB + "0" + TAB + "g1" + TAB + "11" + TAB + "255" + TAB + "3M2N3M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "TTCCTA" + TAB + "``````" + TAB + "AS:i:1" + LS;
  private static final long SAM26_LENGTH = 6 * 6;
  /** SAM File Contents **/
  public static final String SAM26 = "" + SAMHEADER1 + SAM_BODY_26;

  // Similar to CG with internal Ns in cigar - see test 1
  public void test26() throws Exception {
    final String[] args0 = {"-a", "--snps-only"};
    check(REF_SEQS, SAM26, 26, null, 0, SAM26_LENGTH, args0);
  }

  // 12 3456 7890 1234 5678 90
  // aa atcg actg gtca gcta gg
  // cg acag tt
  // cg acag tt
  // cg acag tt
  // cg acag tt
  // g acag ctc
  // g acag ctc
  // g acag ctc
  // g acag ctc
  private static final String SAM27 = "" + SAMHEADER1
      + "1a" + TAB + "0" + TAB + "g1" + TAB + "5" + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACAGTT" + TAB + "````````" + TAB + "AS:i:1" + LS
      + "1b" + TAB + "0" + TAB + "g1" + TAB + "5" + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACAGTT" + TAB + "````````" + TAB + "AS:i:1" + LS
      + "1c" + TAB + "0" + TAB + "g1" + TAB + "5" + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACAGTT" + TAB + "````````" + TAB + "AS:i:1" + LS
      + "1d" + TAB + "0" + TAB + "g1" + TAB + "5" + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACAGTT" + TAB + "````````" + TAB + "AS:i:1" + LS
      + "2a" + TAB + "0" + TAB + "g1" + TAB + "6" + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "GACAGCTC" + TAB + "````````" + TAB + "AS:i:1" + LS
      + "2b" + TAB + "0" + TAB + "g1" + TAB + "6" + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "GACAGCTC" + TAB + "````````" + TAB + "AS:i:1" + LS
      + "2c" + TAB + "0" + TAB + "g1" + TAB + "6" + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "GACAGCTC" + TAB + "````````" + TAB + "AS:i:1" + LS
      + "2d" + TAB + "0" + TAB + "g1" + TAB + "6" + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "GACAGCTC" + TAB + "````````" + TAB + "AS:i:1" + LS;
  private static final long SAM27_LENGTH = 8 * 8;

  // exercise SNP splitting when no indels but complex calling
  public void test27() throws Exception {
    final String[] args0 = {"-a", "--snps-only"};
    check(REF_SEQS, SAM27, 27, null, 0, SAM27_LENGTH, args0);
  }

  // 12 3456 7890 1234 5678 90
  // aa atcg actg gtca gcta gg
  // cg acag tt
  // cg acag tt
  // cg acag tt
  // cg acag tt
  // g acag ttc
  // g acag ttc
  // g acag ttc
  // g acag ttc
  private static final String SAM28 = "" + SAMHEADER1
      + "1a" + TAB + "0" + TAB + "g1" + TAB + "5" + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACAGTT" + TAB + "````````" + TAB + "AS:i:1" + LS
      + "1b" + TAB + "0" + TAB + "g1" + TAB + "5" + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACAGTT" + TAB + "````````" + TAB + "AS:i:1" + LS
      + "1c" + TAB + "0" + TAB + "g1" + TAB + "5" + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACAGTT" + TAB + "````````" + TAB + "AS:i:1" + LS
      + "1d" + TAB + "0" + TAB + "g1" + TAB + "5" + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACAGTT" + TAB + "````````" + TAB + "AS:i:1" + LS
      + "2a" + TAB + "0" + TAB + "g1" + TAB + "6" + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "GACAGTTC" + TAB + "````````" + TAB + "AS:i:1" + LS
      + "2b" + TAB + "0" + TAB + "g1" + TAB + "6" + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "GACAGTTC" + TAB + "````````" + TAB + "AS:i:1" + LS
      + "2c" + TAB + "0" + TAB + "g1" + TAB + "6" + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "GACAGTTC" + TAB + "````````" + TAB + "AS:i:1" + LS
      + "2d" + TAB + "0" + TAB + "g1" + TAB + "6" + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "GACAGTTC" + TAB + "````````" + TAB + "AS:i:1" + LS;
  private static final long SAM28_LENGTH = 8 * 8;

  // exercise SNP splitting when no indels but complex calling
  public void test28() throws Exception {
    final String[] args0 = {"-a", "--snps-only"};
    check(REF_SEQS, SAM28, 28, null, 0, SAM28_LENGTH, args0);
  }

  public void test29() throws Exception {
    final String[] args0 = {"-a", "--snps-only"};
    check(REF_SEQS, SAM27, 29, null, 0, SAM27_LENGTH, args0);
  }

  // exercise SNP splitting when no indels but complex calling
  public void test30() throws Exception {
    final String[] args0 = {"-a", "--snps-only"};
    check(REF_SEQS, SAM28, 30, null, 0, SAM28_LENGTH, args0);
  }

  // exercise SNP splitting when no indels but complex calling
  public void test31() throws Exception {
    final String[] args0 = {"--snps-only"};
    check(REF_SEQS, SAM27, 31, null, 0, SAM27_LENGTH, args0);
  }

  public void test32() throws Exception {
    final String[] args0 = {"--snps-only"};
    check(REF_SEQS, SAM28, 32, null, 0, SAM28_LENGTH, args0);
  }

  // 12 3456 7890 _1234 5678 90
  // aa atcg actg _gtca gcta gg
  // cg actg agt
  // cg actg agt
  // cg actg agt
  // cg actg agt
  // Simple insert - test handling of new all-paths for indels
  private static final String SAM37 = "" + SAMHEADER1
      + "1a" + TAB + "0" + TAB + "g1" + TAB + "5" + TAB + "255" + TAB + "6M1I2M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACTGAGT" + TAB + "`````````" + TAB + "AS:i:1" + LS
      + "1b" + TAB + "0" + TAB + "g1" + TAB + "5" + TAB + "255" + TAB + "6M1I2M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACTGAGT" + TAB + "`````````" + TAB + "AS:i:1" + LS
      + "1c" + TAB + "0" + TAB + "g1" + TAB + "5" + TAB + "255" + TAB + "6M1I2M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACTGAGT" + TAB + "`````````" + TAB + "AS:i:1" + LS
      + "1d" + TAB + "0" + TAB + "g1" + TAB + "5" + TAB + "255" + TAB + "6M1I2M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACTGAGT" + TAB + "`````````" + TAB + "AS:i:1" + LS
      + "1e" + TAB + "0" + TAB + "g1" + TAB + "5" + TAB + "255" + TAB + "6M1I2M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACTGAGT" + TAB + "`````````" + TAB + "AS:i:1" + LS;
  private static final long SAM37_LENGTH = 5 * 9;

  // exercise SNP splitting when no indels but complex calling
  public void test37() throws Exception {
    final String[] args0 = {"-a", "--Xno-complex-calls"};

    //remove g1 11  x i x 0.0 5 0.050 A 5 0.050 for singleton caller
    check(REF_SEQS, SAM37, 37, null, 0, SAM37_LENGTH, args0);
  }

  public void test38() throws Exception {
    final String[] args0 = {"-a"};
    check(REF_SEQS, SAM37, 38, null, 0, SAM37_LENGTH, args0);
  }

  // private static final String SAM41 = ""
  // + SAMHEADER1
  // + "1a" + TAB + "0" + TAB + "g1" + TAB + "5" + TAB + "255" + TAB + "6M1D2M"
  // + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACTGTC" + TAB + "````````" +
  // TAB + "AS:i:1" + LS
  // + "1b" + TAB + "0" + TAB + "g1" + TAB + "5" + TAB + "255" + TAB + "6M1D2M"
  // + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACTGTC" + TAB + "````````" +
  // TAB + "AS:i:1" + LS
  // + "1c" + TAB + "0" + TAB + "g1" + TAB + "5" + TAB + "255" + TAB + "6M1D2M"
  // + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACTGTC" + TAB + "````````" +
  // TAB + "AS:i:1" + LS
  // + "1d" + TAB + "0" + TAB + "g1" + TAB + "5" + TAB + "255" + TAB + "6M1D2M"
  // + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACTGTC" + TAB + "````````" +
  // TAB + "AS:i:1" + LS
  // + "1e" + TAB + "0" + TAB + "g1" + TAB + "5" + TAB + "255" + TAB + "6M1D2M"
  // + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACTGTC" + TAB + "````````" +
  // TAB + "AS:i:1" + LS
  // ;
  // public void test41() {
  // final String[] args0 = {
  // "-a", "--indels", "--Xcomplex-indels"//, "--complex-calls"
  // };
  // check(REF_SEQS, SAM41, args0, 41, "", 0);
  // }

  // 12 3456 789_0 1_234 5678 90
  // aa atcg act_g g_tca gcta gg
  // cg actgg g_t
  // cg actgg g_t
  // cg actgg g_t
  // cg act_g ggt
  // cg act_g ggt
  // cg act_g ggt
  // Simple insert - right at the end of an interesting region
  // private static final String SAM39 = ""
  // + SAMHEADER1
  // + "1a" + TAB + "0" + TAB + "g1" + TAB + "5" + TAB + "255" + TAB + "5=1I3="
  // + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACTGGGT" + TAB + "`````````"
  // + TAB + "AS:i:1" + LS
  // + "1b" + TAB + "0" + TAB + "g1" + TAB + "5" + TAB + "255" + TAB + "5=1I3="
  // + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACTGGGT" + TAB + "`````````"
  // + TAB + "AS:i:1" + LS
  // + "1c" + TAB + "0" + TAB + "g1" + TAB + "5" + TAB + "255" + TAB + "5=1I3="
  // + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACTGGGT" + TAB + "`````````"
  // + TAB + "AS:i:1" + LS
  // + "2a" + TAB + "0" + TAB + "g1" + TAB + "5" + TAB + "255" + TAB + "7=1I1="
  // + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACTGGGT" + TAB + "`````````"
  // + TAB + "AS:i:1" + LS
  // + "2b" + TAB + "0" + TAB + "g1" + TAB + "5" + TAB + "255" + TAB + "7=1I1="
  // + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACTGGGT" + TAB + "`````````"
  // + TAB + "AS:i:1" + LS
  // + "2c" + TAB + "0" + TAB + "g1" + TAB + "5" + TAB + "255" + TAB + "7=1I1="
  // + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACTGGGT" + TAB + "`````````"
  // + TAB + "AS:i:1" + LS
  // ;
  // exercise SNP splitting when no indels but complex calling
  // public void test39() {
  // final String[] args0 = {
  // "-a", "--indels", "--complex-calls"
  // };
  // check(REF_SEQS, SAM39, args0, 39, "", 0);
  // }

  /** Reference Sequence **/
  public static final String REF_SEQS_MULTI = "" + ">g1" + LS + "aa" + "atcg" + "acac" + "gtca" + "gcta" + "gg" + LS + ">g2" + LS + "aa" + "atcg" + "acac" + "gtca" + "gcta" + "gg" + LS;
  /*
   * 123456789012345 REF aaatcga*actgacgt atcgacac atcgacac atcgacac
   */
  static final String SAM_TEST_MULTI_SEQUENCE = ""
      + "@HD" + TAB + "VN:1.0" + TAB + "SO:coordinate" + LS
      + "@SQ" + TAB + "SN:g1" + TAB + "LN:20" + LS
      + "@SQ" + TAB + "SN:g2" + TAB + "LN:20" + LS
      + TEMPLATE_SDF_ID + LS
      + "a0" + TAB + "0" + TAB + "g1" + TAB + "3" + TAB + "255" + TAB + "4=1X1I2=" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ATCGGCCA" + TAB + "````````" + TAB + "AS:i:0" + LS
      + "a1" + TAB + "0" + TAB + "g1" + TAB + "3" + TAB + "255" + TAB + "4=1X1I2=" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ATCGGCCA" + TAB + "````````" + TAB + "AS:i:0" + LS
      + "a2" + TAB + "0" + TAB + "g1" + TAB + "3" + TAB + "255" + TAB + "4=1X1I2=" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ATCGGCCA" + TAB + "````````" + TAB + "AS:i:0" + LS
      + "a3" + TAB + "0" + TAB + "g1" + TAB + "3" + TAB + "255" + TAB + "4=1X1I2=" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ATCGGCCA" + TAB + "````````" + TAB + "AS:i:0" + LS
      // Added the following line to test the ignoring of zero-length reads rather than creating error messages from picard
      + "a4" + TAB + "4" + TAB + "g1" + TAB + "3" + TAB + "0" + TAB + "*" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "*" + TAB + "*" + LS
      + "b0" + TAB + "0" + TAB + "g2" + TAB + "3" + TAB + "255" + TAB + "4=1X1I2=" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ATCGTGCA" + TAB + "````````" + TAB + "AS:i:0" + LS
      + "b1" + TAB + "0" + TAB + "g2" + TAB + "3" + TAB + "255" + TAB + "4=1X1I2=" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ATCGTGCA" + TAB + "````````" + TAB + "AS:i:0" + LS
      + "b2" + TAB + "0" + TAB + "g2" + TAB + "3" + TAB + "255" + TAB + "4=1X1I2=" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ATCGTGCA" + TAB + "````````" + TAB + "AS:i:0" + LS
      + "b3" + TAB + "0" + TAB + "g2" + TAB + "3" + TAB + "255" + TAB + "4=1X1I2=" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ATCGTGCA" + TAB + "````````" + TAB + "AS:i:0" + LS;
  private static final long SAM_TMS_LENGTH = 8 * 8;

  public void test40() throws Exception {
    check(REF_SEQS_MULTI, SAM_TEST_MULTI_SEQUENCE, 40, null, 0, SAM_TMS_LENGTH);
  }

  /** Reference Sequence **/
  public static final String COMPLEX_EQUIVALENT_TEMPLATE = "" + ">g1" + LS + "acccgtagagagagatgctatt" + LS;
  /*
   * 12345678901__2345678 REF acccgt__agagagaga__tgctatt cgt__agagagagagatgct
   * cgt__agagagagagatgct cgt__agagagagagatgct cgt__agagagagagatgct
   * cgtagagagagaga__tgct cgtagagagagaga__tgct cgtagagagagaga__tgct
   * cgtagagagagaga__tgct
   */
  static final String COMPLEX_EQUIV_READ = "cgtagagagagagatgct";
  private static final String COMPLEX_EQUIV_QUALITY = "``````````````````";

  static final String SAM_COMPLEX_EQUIV_HEADER = ""
      + "@HD" + TAB + "VN:1.0" + TAB + "SO:coordinate" + LS
      + "@SQ" + TAB + "SN:g1" + TAB + "LN:22" + LS
      + TEMPLATE_SDF_ID + LS;
  static final String SAM_COMPLEX_EQUIV_BASE_LINE = "%s\t0\tg1\t4\t255\t%s\t*\t0\t0\t%s\t%s\tAS:i:3%n";
  static final String SAM_TEST_COMPLEX_EQUIVALENT = ""
    + SAM_COMPLEX_EQUIV_HEADER
    + String.format(SAM_COMPLEX_EQUIV_BASE_LINE, "a0", "3=2I13=", COMPLEX_EQUIV_READ, COMPLEX_EQUIV_QUALITY)
    + String.format(SAM_COMPLEX_EQUIV_BASE_LINE, "a1", "3=2I13=", COMPLEX_EQUIV_READ, COMPLEX_EQUIV_QUALITY)
    + String.format(SAM_COMPLEX_EQUIV_BASE_LINE, "a2", "3=2I13=", COMPLEX_EQUIV_READ, COMPLEX_EQUIV_QUALITY)
    + String.format(SAM_COMPLEX_EQUIV_BASE_LINE, "a3", "3=2I13=", COMPLEX_EQUIV_READ, COMPLEX_EQUIV_QUALITY)
    + String.format(SAM_COMPLEX_EQUIV_BASE_LINE, "b0", "12=2I4=", COMPLEX_EQUIV_READ, COMPLEX_EQUIV_QUALITY)
    + String.format(SAM_COMPLEX_EQUIV_BASE_LINE, "b1", "12=2I4=", COMPLEX_EQUIV_READ, COMPLEX_EQUIV_QUALITY)
    + String.format(SAM_COMPLEX_EQUIV_BASE_LINE, "b2", "12=2I4=", COMPLEX_EQUIV_READ, COMPLEX_EQUIV_QUALITY)
    + String.format(SAM_COMPLEX_EQUIV_BASE_LINE, "b3", "12=2I4=", COMPLEX_EQUIV_READ, COMPLEX_EQUIV_QUALITY);

  private static final long SAM_TCE_LENGTH = 8 * 18;

  public void test41() throws Exception {
    check(COMPLEX_EQUIVALENT_TEMPLATE, SAM_TEST_COMPLEX_EQUIVALENT, 41, null, 0, SAM_TCE_LENGTH);
  }
//                                                  "cgtagagagagagatgct"
//                                               "acccgtagagagagatgctatt"
  static final String COMPLEX_EQUIV_READ_EQUALS_A = "===ag=============";
  static final String COMPLEX_EQUIV_READ_EQUALS_B = "============ga====";

  static final String SAM_TEST_COMPLEX_EQUIVALENT_EQUALS = ""
    + SAM_COMPLEX_EQUIV_HEADER
    + String.format(SAM_COMPLEX_EQUIV_BASE_LINE, "a0", "3=2I13=", COMPLEX_EQUIV_READ_EQUALS_A, COMPLEX_EQUIV_QUALITY)
    + String.format(SAM_COMPLEX_EQUIV_BASE_LINE, "a1", "3=2I13=", COMPLEX_EQUIV_READ_EQUALS_A, COMPLEX_EQUIV_QUALITY)
    + String.format(SAM_COMPLEX_EQUIV_BASE_LINE, "a2", "3=2I13=", COMPLEX_EQUIV_READ_EQUALS_A, COMPLEX_EQUIV_QUALITY)
    + String.format(SAM_COMPLEX_EQUIV_BASE_LINE, "a3", "3=2I13=", COMPLEX_EQUIV_READ_EQUALS_A, COMPLEX_EQUIV_QUALITY)
    + String.format(SAM_COMPLEX_EQUIV_BASE_LINE, "b0", "12=2I4=", COMPLEX_EQUIV_READ_EQUALS_B, COMPLEX_EQUIV_QUALITY)
    + String.format(SAM_COMPLEX_EQUIV_BASE_LINE, "b1", "12=2I4=", COMPLEX_EQUIV_READ_EQUALS_B, COMPLEX_EQUIV_QUALITY)
    + String.format(SAM_COMPLEX_EQUIV_BASE_LINE, "b2", "12=2I4=", COMPLEX_EQUIV_READ_EQUALS_B, COMPLEX_EQUIV_QUALITY)
    + String.format(SAM_COMPLEX_EQUIV_BASE_LINE, "b3", "12=2I4=", COMPLEX_EQUIV_READ_EQUALS_B, COMPLEX_EQUIV_QUALITY);

  public void test41SamEquals() throws Exception {
    check(COMPLEX_EQUIVALENT_TEMPLATE, SAM_TEST_COMPLEX_EQUIVALENT_EQUALS, 41, null, 0, SAM_TCE_LENGTH);
  }

  /*
   * 12345678901__2345678 REF acccgt__agagagaga__tgctatt cgtagagagcgaga__tgct
   * cgtagagagcgaga__tgct cgtagagagcgaga__tgct cgtagagagcgaga__tgct
   * cgt__agagcgagagatgct cgt__agagcgagagatgct cgt__agagcgagagatgct
   * cgt__agagcgagagatgct
   */
  static final String SAM_TEST_COMPLEX_EQUIVALENT_2 = ""
      + "@HD" + TAB + "VN:1.0" + TAB + "SO:coordinate" + LS
      + TEMPLATE_SDF_ID + LS
      + "@SQ" + TAB + "SN:g1" + TAB + "LN:22" + LS
      + "a0" + TAB + "0" + TAB + "g1" + TAB + "4" + TAB + "255" + TAB + "3=2I4=1X8=" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "cgtagagagcgagatgct" + TAB + "``````````````````" + TAB + "AS:i:3" + LS + "a1" + TAB + "0" + TAB + "g1" + TAB + "4" + TAB + "255" + TAB + "3=2I4=1X8=" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "cgtagagagcgagatgct" + TAB + "``````````````````" + TAB + "AS:i:3" + LS + "a2" + TAB + "0" + TAB + "g1" + TAB + "4" + TAB + "255" + TAB + "3=2I4=1X8=" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "cgtagagagcgagatgct" + TAB + "``````````````````" + TAB + "AS:i:3" + LS + "a3" + TAB + "0" + TAB + "g1" + TAB + "4" + TAB + "255" + TAB + "3=2I4=1X8=" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "cgtagagagcgagatgct" + TAB + "``````````````````" + TAB + "AS:i:3" + LS + "b0" + TAB + "0" + TAB + "g1" + TAB + "4" + TAB + "255" + TAB + "7=1X4=2I4=" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "cgtagagcgagagatgct" + TAB + "``````````````````" + TAB + "AS:i:3" + LS + "b1" + TAB + "0" + TAB + "g1" + TAB + "4" + TAB + "255" + TAB + "7=1X4=2I4=" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "cgtagagcgagagatgct" + TAB + "``````````````````" + TAB + "AS:i:3" + LS + "b2" + TAB + "0" + TAB + "g1" + TAB + "4" + TAB + "255" + TAB + "7=1X4=2I4=" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "cgtagagcgagagatgct" + TAB + "``````````````````" + TAB + "AS:i:3" + LS + "b3" + TAB + "0" + TAB + "g1" + TAB + "4" + TAB + "255" + TAB + "7=1X4=2I4=" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "cgtagagcgagagatgct" + TAB + "``````````````````" + TAB + "AS:i:3" + LS;

  public void testComplexEquivalentWithBrokenSymmetry() throws Exception {
    final String[] args0 = {};
    check(COMPLEX_EQUIVALENT_TEMPLATE, SAM_TEST_COMPLEX_EQUIVALENT_2, 42, null, 0, 144L/*regression*/, args0);
  }

  // 12 3456 7890 1234 5678 90
  // aa atcg actg gtca gcta gg
  // cg accg gt
  // cg accg gt
  // cg accg gt
  // cg accg gt
  // Simple snp in an ambiguous region
  private static final String SAM_AMBIGUOUS = "" + SAMHEADER1 + "1a" + TAB + "0" + TAB + "g1" + TAB + "5" + TAB + "255" + TAB + "4=1X3=" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACCGGT" + TAB + "````````" + TAB + "AS:i:1" + LS + "1b" + TAB + "0" + TAB + "g1" + TAB + "5" + TAB + "255" + TAB + "4=1X3=" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACCGGT" + TAB + "````````" + TAB + "AS:i:1" + LS + "1c" + TAB + "0" + TAB + "g1" + TAB + "5" + TAB + "255" + TAB + "4=1X3=" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACCGGT" + TAB + "````````" + TAB + "AS:i:1" + LS + "1d" + TAB + "0" + TAB + "g1" + TAB + "5" + TAB + "255" + TAB + "4=1X3=" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACCGGT" + TAB + "````````" + TAB + "AS:i:1" + LS + "1e" + TAB + "0" + TAB + "g1" + TAB + "5" + TAB + "255" + TAB + "4=1X3=" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACCGGT" + TAB + "````````" + TAB + "AS:i:1" + TAB + "IH:i:2" + LS;

  public void test43AmbiguityFilter() throws Exception {
    final String[] args0 = {"--filter-ambiguity", "10%"};
    check(REF_SEQS, SAM_AMBIGUOUS, 43, null, 0, 40L/*regression*/, args0);
  }

  // 12 3456 7890 1234 5678 90
  // aa atcg actg gtca gcta gg
  // atcg accg
  // cg accg gt
  // cg accg gt
  // cg accg gt
  // cg accg gt
  // accg gtca
  // Test handling of coverage threshold
  private static final String SAM_COMPLEX_BODY = ""
      + "0a" + TAB + "0" + TAB + "g1" + TAB + "3" + TAB + "255" + TAB + "5=2X1=" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ATCGATCG" + TAB + "````````" + TAB + "AS:i:1" + LS
      + "1a" + TAB + "0" + TAB + "g1" + TAB + "5" + TAB + "255" + TAB + "3=2X3=" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGATCGGT" + TAB + "````````" + TAB + "AS:i:1" + LS
      + "1b" + TAB + "0" + TAB + "g1" + TAB + "5" + TAB + "255" + TAB + "3=2X3=" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGATCGGT" + TAB + "````````" + TAB + "AS:i:1" + LS
      + "1c" + TAB + "0" + TAB + "g1" + TAB + "5" + TAB + "255" + TAB + "3=2X3=" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGATCGGT" + TAB + "````````" + TAB + "AS:i:1" + LS
      + "1d" + TAB + "0" + TAB + "g1" + TAB + "5" + TAB + "255" + TAB + "3=2X3=" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGATCGGT" + TAB + "````````" + TAB + "AS:i:1" + LS
      + "1e" + TAB + "0" + TAB + "g1" + TAB + "5" + TAB + "255" + TAB + "3=2X3=" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGATCGGT" + TAB + "````````" + TAB + "AS:i:1" + LS
      + "2a" + TAB + "0" + TAB + "g1" + TAB + "7" + TAB + "255" + TAB + "1=2X5=" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ATCGGTCA" + TAB + "````````" + TAB + "AS:i:1" + LS;
  private static final String SAM_COMPLEX_COVERAGE = "" + SAMHEADER1 + SAM_COMPLEX_BODY;

  public void test44ComplexCoverageThreshold() throws Exception {
    final String[] args0 = {"--filter-depth", "6", "--all"};
    check(REF_SEQS, SAM_COMPLEX_COVERAGE, 44, null, 0, 56L/*regression*/, args0);
  }

  /** SAM format file after mapping and alignment **/
  public static final String SAM_BODY_45 = "" + "0" + TAB + "0" + TAB + "g1" + TAB + "3" + TAB + "5" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ATCGACTG" + TAB + "````````" + TAB + "AS:i:0" + TAB + "RG:Z:RG1" + LS + "1" + TAB + "0" + TAB + "g1" + TAB + "3" + TAB + "5" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ATCGACTG" + TAB + "````````" + TAB + "AS:i:0" + TAB + "RG:Z:RG1" + LS + "2" + TAB + "0" + TAB + "g1" + TAB + "5" + TAB + "5" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACTGTT" + TAB + "````````" + TAB + "AS:i:1" + TAB + "RG:Z:RG1" + LS + "3" + TAB + "0" + TAB + "g1" + TAB + "6" + TAB + "5" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "GACTGCTC" + TAB + "````````" + TAB + "AS:i:1" + TAB + "RG:Z:RG1" + LS + "4" + TAB + "0" + TAB + "g1" + TAB + "6" + TAB + "5" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "GACTGCTC" + TAB + "````````" + TAB + "AS:i:1" + TAB + "RG:Z:RG1" + LS + "5" + TAB + "0" + TAB + "g1" + TAB + "11" + TAB + "5" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "TTCAGCTA" + TAB + "````````" + TAB + "AS:i:1" + TAB + "RG:Z:RG1" + LS;

  /** SAM File Contents **/
  public static final String SAM_45 = "" + SAMHEADER1 + SAM_BODY_45;

  // test that new mapq scores being actually used and reflected in statistics
  // and posterior - see test1
  public void test45() throws Exception {
    final String[] args0 = {"-a", "--Xno-complex-calls", "--snps-only"};
    check(REF_SEQS, SAM_45, 45, null, 0, 48L/*regression*/, args0);
  }

  /** SAM format file after mapping and alignment **/
  public static final String SAM_BODY_46 = "" + "0" + TAB + "0" + TAB + "g1" + TAB + "3" + TAB + "55" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ATCGACTG" + TAB + "````````" + TAB + "AS:i:0" + TAB + "RG:Z:RG1" + LS + "1" + TAB + "0" + TAB + "g1" + TAB + "3" + TAB + "55" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ATCGACTG" + TAB + "````````" + TAB + "AS:i:0" + TAB + "RG:Z:RG1" + LS + "2" + TAB + "0" + TAB + "g1" + TAB + "5" + TAB + "55" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACTGTT" + TAB + "````````" + TAB + "AS:i:1" + TAB + "RG:Z:RG1" + LS + "3" + TAB + "0" + TAB + "g1" + TAB + "6" + TAB + "55" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "GACTGCTC" + TAB + "````````" + TAB + "AS:i:1" + TAB + "RG:Z:RG1" + LS + "4" + TAB + "0" + TAB + "g1" + TAB + "6" + TAB + "55" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "GACTGCTC" + TAB + "````````" + TAB + "AS:i:1" + TAB + "RG:Z:RG1" + LS + "5" + TAB + "0" + TAB + "g1" + TAB + "11" + TAB + "55" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "TTCAGCTA" + TAB + "````````" + TAB + "AS:i:1" + TAB + "RG:Z:RG1" + LS;

  /** SAM File Contents **/
  public static final String SAM_46 = "" + SAMHEADER1 + SAM_BODY_46;

  // test that new mapq scores being actually used and reflected in statistics
  // and posterior - see test1
  public void test46() throws Exception {
    final String[] args0 = {"-a", "--Xno-complex-calls", "--snps-only"};
    check(REF_SEQS, SAM_46, 46, null, 0, 48L/*regression*/, args0);
  }

  // test output 'a' for ambiguous - see test1
  public void test47() throws Exception {
    final String[] args0 = {"-a", "--Xno-complex-calls", "--snps-only", "--filter-ambiguity", "10%"};
    check(REF_SEQS, SAM1_AMB, 47, null, 0, 48L/*regression*/, args0);
  }

  // test output 'c' for over coverage threshold - see test1
  public void test48() throws Exception {
    final String[] args0 = {"-a", "--Xno-complex-calls", "--snps-only", "--filter-depth=3"};
    check(REF_SEQS, SAM1, 48, null, 0, 48L/*regression*/, args0);
  }

  private static final String SAM_HEADER = "@HD\tVN:1.4\tSO:coordinate\n";

  private void checkNoDefaultMachineErrors(final String refSeq, final String sam, final String[] args0, final String errorMsg, final String logMsg, final int errCode) throws Exception {
    try (final TestDirectory dir = new TestDirectory("testcheck")) {
      final File output = new File(dir, "variant_out");
      final File input = new File(dir, "variant_in");
      if (!input.mkdir()) {
        throw new IOException("Whut, couldn't make test directory");
      }
      final File insam = new File(input, OUT_SAM + ".gz");
      BgzipFileHelper.streamToBgzipFile(new ByteArrayInputStream((SAM_HEADER + sam).getBytes()), insam);
      new TabixIndexer(insam, new File(input, OUT_SAM + ".gz.tbi")).saveSamIndex();
      final String outn = output.getPath();
      final File templ = ReaderTestUtils.getDNADir(refSeq);
      try {
        final String[] args = Utils.append(args0, "-t", templ.getPath(), "-o", outn, "--Xpriors", "testhumanprior", "-Z", insam.getPath(), "--" + AbstractMultisampleCli.NO_CALIBRATION);
        final MainResult r = MainResult.run(new SingletonCli(), args);
        assertEquals(r.err(), errCode, r.rc());
        assertTrue(r.err(), r.err().contains(errorMsg));
        if (logMsg != null && !logMsg.equals("")) {
          final String lout = FileUtils.fileToString(new File(outn, "snp.log"));
          assertTrue(lout.contains(logMsg));
        }
      } finally {
        FileHelper.deleteAll(templ);
      }
    }
  }

  public void testMachineErrorDetection() throws Exception {
    checkNoDefaultMachineErrors(">a\nacgta", "", new String[] {}, "No read groups found. Unable to determine machine error rate.", null, 1);
    checkNoDefaultMachineErrors(">a\nacgta", "@SQ\tSN:a\tLN:5\n@RG\tSM:foo\tID:sounique\tPL:ILLUMINA\n", new String[] {}, "", "Loading machine errors for: illumina", 0);
    checkNoDefaultMachineErrors(">a\nacgta", "@SQ\tSN:a\tLN:5\n@RG\tSM:foo\tID:sounique\tPL:COMPLETE\n", new String[] {}, "", "Loading machine errors for: complete", 0);
    checkNoDefaultMachineErrors(">a\nacgta", "@SQ\tSN:a\tLN:5\n@RG\tSM:foo\tID:sounique\tPL:LS454\n", new String[] {}, "", "Loading machine errors for: ls454_se", 0);
    checkNoDefaultMachineErrors(">a\nacgta", "@SQ\tSN:a\tLN:5\n@RG\tSM:foo\tID:sounique\tPL:LS454\tPI:500\n", new String[] {}, "", "Loading machine errors for: ls454_pe", 0);
  }

  public void testAllSequenceRestriction() throws Exception {
    checkReaderRestriction("g2", "variantnano_tabix2.vcf");
  }

  //region handling
  public void testAllReaderRestriction() throws Exception {
    checkReaderRestriction("g1:9-18", "variantnano_tabix.vcf");
  }

  public void checkReaderRestriction(String restriction, String results) throws Exception {
    // requires a tabix file to be able to do the reader restriction
    // check(REF_SEQS11, SAM11, args0, 1, "", 0);
    try (final TestDirectory dir = new TestDirectory("testcheck")) {
      final File output = new File(dir, "variant_out");
      final File input = new File(dir, "variant_in");
      if (!input.mkdir()) {
        throw new IOException("Whut, couldn't make test directory");
      }
      final File sam = new File(input, "sam11.sam.gz");
      FileHelper.resourceToFile("com/rtg/variant/resources/sam11.sam.gz", sam);
      //System.err.println(FileHelper.gzFileToString(sam));
      FileHelper.resourceToFile("com/rtg/variant/resources/sam11.sam.gz.tbi", new File(input, "sam11.sam.gz.tbi"));
      final String outn = output.getPath();
      final String inn = sam.getPath();
      final File templ = ReaderTestUtils.getDNADir(SharedSamConstants.REF_SEQS11);
      try {
        final MainResult r = MainResult.run(new SingletonCli(), "-t", templ.getPath(), "-o", outn, "--Xpriors", "testhumanprior", "-Z", inn, "-a", "--region", restriction, "--snps-only", "-m", "default", "--keep-duplicates", "--" + AbstractMultisampleCli.NO_CALIBRATION);
        assertEquals(r.err(), 0, r.rc());
        final String result = FileUtils.fileToString(new File(output, VariantParams.VCF_OUT_SUFFIX));
        final String actualFixed = TestUtils.sanitizeVcfHeader(result);

        //assertEquals(FileHelper.resourceToString("com/rtg/variant/resources/" + results), actualFixed);
        mNano.check(results, actualFixed, false);
      } finally {
        FileHelper.deleteAll(templ);
      }
    }
  }

  public void testHyperComplexRegion() throws Exception {
    try (final TestDirectory dir = new TestDirectory("test")) {
      final File template = ReaderTestUtils.getDNADir(FileHelper.resourceToString("com/rtg/variant/resources/hypertemplate.fa"), new File(dir, "template"));
      final File sam = new File(dir, "sam.sam.gz");
      VariantTestUtils.bgzipAndIndex(FileHelper.resourceToString("com/rtg/variant/resources/hypersam.sam"), sam);
      final File output = new File(dir, "variant_out");
      final MainResult r = MainResult.run(new SingletonCli(), "-t", template.getPath(), "-o", output.getPath(), "-m", "illumina", "--keep-duplicates", sam.getPath(), "--" + AbstractMultisampleCli.NO_CALIBRATION);
      assertEquals(r.err(), 0, r.rc());
      final String result = StringUtils.grep(FileHelper.gzFileToString(new File(output, "snps.vcf.gz")), "^[^#]").replaceAll("\n|\r\n", LS);

      //turn on following in singleton caller case
      mNano.check("hyper_oldcli.vcf", result, false);
    }
  }

  static final String SAMHEADER11 = ""
       + "@HD" + TAB + "VN:1.0" + TAB + "SO:coordinate" + LS
       + "@SQ" + TAB + "SN:g1" + TAB + "LN:20" + LS
       + "@SQ" + TAB + "SN:gempty" + TAB + "LN:0" + LS
       + TEMPLATE_SDF_ID + LS
       + "@RG" + TAB + "ID:RG1" + TAB + "SM:TEST" + TAB + "PL:ILLUMINA" + LS;

  // 12 3456 7890 1234 5678 90
  // aa atcg actg gtca gcta gg
  // cg actg tca
  // cg actg tca
  // cg actg tca
  // cg actg tca
  // Simple delete
  private static final String SAM49 = ""
      + SAMHEADER11
      + "1a" + TAB + "0" + TAB + "g1" + TAB + "5" + TAB + "255" + TAB + "6M1D3M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACTGTCA" + TAB + "`````````" + TAB + "AS:i:1" + LS
      + "1b" + TAB + "0" + TAB + "g1" + TAB + "5" + TAB + "255" + TAB + "6M1D3M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACTGTCA" + TAB + "`````````" + TAB + "AS:i:1" + LS
      + "1c" + TAB + "0" + TAB + "g1" + TAB + "5" + TAB + "255" + TAB + "6M1D3M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACTGTCA" + TAB + "`````````" + TAB + "AS:i:1" + LS
      + "1d" + TAB + "0" + TAB + "g1" + TAB + "5" + TAB + "255" + TAB + "6M1D3M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACTGTCA" + TAB + "`````````" + TAB + "AS:i:1" + LS
      + "1e" + TAB + "0" + TAB + "g1" + TAB + "5" + TAB + "255" + TAB + "6M1D3M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACTGTCA" + TAB + "`````````" + TAB + "AS:i:1" + LS;
  private static final int SAM49_LENGTH = 5 * 9;

  public void test49() throws Exception {
    check(REF_SEQS, SAM49, 49, null, 0, SAM49_LENGTH, "-a");
  }

  public void test50() throws Exception {
    check(REF_SEQS, SAM49, 50, null, "", 0, 10, null, true, SAM49_LENGTH, "-a");
  }

  public void test51() throws Exception {
    final String sam = FileHelper.resourceToString("com/rtg/variant/resources/vcftest.sam");
    final String tmplfa = FileHelper.resourceToString("com/rtg/variant/resources/vcftesttmpl.fa");
    check(tmplfa, sam, 51, null, "", 0, 2, null, true, 172L/*regression*/);
  }

  public void test52() throws Exception {
    check(REF_SEQS, SAM1_AMB, 52, null, "", 0, 10, null, true, 48L/*regression*/, "-a", "--snps-only", "--filter-ambiguity", "10%");
  }

  public void test53() throws Exception { // same as test17 except haploid
    check(REF_SEQS, SAM1_I2, 53, null, "Processing g1", 0, 1, SharedSamConstants.REF_HAPLOID, true, SAM1_I2_LENGTH, "-a", "--sex=either");
  }

  public void test53PloidyOverride() throws Exception { // same as test17 and test53 but via ploidy override
    check(REF_SEQS, SAM1_I2, 17, null, "Processing g1", 0, 1, SharedSamConstants.REF_HAPLOID, true, SAM1_I2_LENGTH, "-a", "--ploidy=diploid");
    check(REF_SEQS, SAM1_I2, 53, null, "Processing g1", 0, 1, SharedSamConstants.REF_DIPLOID, true, SAM1_I2_LENGTH, "-a", "--ploidy=haploid");
  }

  //similar to test6 but long enough to get multiple chunks
  private static final String SAM54 = ""
      + "@HD" + TAB + "VN:1.0" + TAB + "SO:coordinate" + LS + "@SQ" + TAB + "SN:s1" + TAB + "LN:5000" + LS
      + TEMPLATE_SDF_ID + LS
      + "0"    + TAB + "0" + TAB + "s1" + TAB + "3"    + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ATCGACTG" + TAB + "*" + TAB + "AS:i:0" + LS
      + "1"    + TAB + "0" + TAB + "s1" + TAB + "3"    + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ATCGACTG" + TAB + "*" + TAB + "AS:i:0" + LS
      + "2"    + TAB + "0" + TAB + "s1" + TAB + "5"    + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACTGTT" + TAB + "*" + TAB + "AS:i:1" + LS
      + "3"    + TAB + "0" + TAB + "s1" + TAB + "6"    + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "GACTGCTC" + TAB + "*" + TAB + "AS:i:1" + LS
      + "4"    + TAB + "0" + TAB + "s1" + TAB + "6"    + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "GACTGCTC" + TAB + "*" + TAB + "AS:i:1" + LS
      + "5"    + TAB + "0" + TAB + "s1" + TAB + "11"   + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "TTCAGCTA" + TAB + "*" + TAB + "AS:i:1" + LS
      + "1000" + TAB + "0" + TAB + "s1" + TAB + "1003" + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ATCGACTG" + TAB + "*" + TAB + "AS:i:0" + LS
      + "1001" + TAB + "0" + TAB + "s1" + TAB + "1003" + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ATCGACTG" + TAB + "*" + TAB + "AS:i:0" + LS
      + "1002" + TAB + "0" + TAB + "s1" + TAB + "1005" + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACTGTT" + TAB + "*" + TAB + "AS:i:1" + LS
      + "1003" + TAB + "0" + TAB + "s1" + TAB + "1006" + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "GACTGCTC" + TAB + "*" + TAB + "AS:i:1" + LS
      + "1004" + TAB + "0" + TAB + "s1" + TAB + "1006" + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "GACTGCTC" + TAB + "*" + TAB + "AS:i:1" + LS
      + "1005" + TAB + "0" + TAB + "s1" + TAB + "1011" + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "TTCAGCTA" + TAB + "*" + TAB + "AS:i:1" + LS
      + "4000" + TAB + "0" + TAB + "s1" + TAB + "4003" + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ATCGACTG" + TAB + "*" + TAB + "AS:i:0" + LS
      + "4001" + TAB + "0" + TAB + "s1" + TAB + "4003" + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ATCGACTG" + TAB + "*" + TAB + "AS:i:0" + LS
      + "4002" + TAB + "0" + TAB + "s1" + TAB + "4005" + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACTGTT" + TAB + "*" + TAB + "AS:i:1" + LS
      + "4003" + TAB + "0" + TAB + "s1" + TAB + "4006" + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "GACTGCTC" + TAB + "*" + TAB + "AS:i:1" + LS
      + "4004" + TAB + "0" + TAB + "s1" + TAB + "4006" + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "GACTGCTC" + TAB + "*" + TAB + "AS:i:1" + LS
      + "4005" + TAB + "0" + TAB + "s1" + TAB + "4011" + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "TTCAGCTA" + TAB + "*" + TAB + "AS:i:1" + LS
      ;
  private static final long SAM54_LENGTH = 18 * 8;

  /** Reference Sequence **/
  public static final String REF_SEQS54;
  static {
    final StringBuilder sb = new StringBuilder();
    sb.append(">s1").append(LS);
    while (sb.length() < 5000) {
      sb.append("aa" + "atcg" + "actg" + "gtca" + "gcta" + "gg");
    }
    sb.append(LS);
    REF_SEQS54 = sb.toString();
  }

  public void test54() throws Exception {
    check(REF_SEQS54, SAM54, 54, null, 0, SAM54_LENGTH);
  }

  // bug 1332
  public void testIdentityInsert() throws Exception {
    try (final TestDirectory dir = new TestDirectory("test")) {
      final File template = ReaderTestUtils.getDNADir(FileHelper.resourceToString("com/rtg/variant/resources/ident-insert2-genome.fa"), new File(dir, "template"));
      final File sam = VariantTestUtils.bgzipAndIndexResource("com/rtg/variant/resources/ident-insert2.sam", dir);
      //      final String samString = FileHelper.resourceToString("com/rtg/variant/resources/ident-insert2.sam");
      //      final File sam = new File(dir, "sam.sam.gz");
      //      bgzipAndIndex(samString, sam);
      final File output = new File(dir, "variant_out");

      final MainResult r = MainResult.run(new SingletonCli(), "-t", template.getPath(),
        "-o", output.getPath(), sam.getPath(),
        "--keep-duplicates", "--" + AbstractMultisampleCli.NO_CALIBRATION);
      assertEquals(r.err(), 0, r.rc());

      final String result = StringUtils.grep(FileHelper.gzFileToString(new File(output, "snps.vcf.gz")), "^[^#]").trim();
      mNano.check("variantnanotest-identityinsert.vcf", result, false);
      //assertTrue(result, result.isEmpty());
    }
  }

  public void testComplexDeletesInHyperComplex() throws Exception {
    try (final TestDirectory dir = new TestDirectory("test")) {
      final File template = ReaderTestUtils.getDNADir(FileHelper.resourceToString("com/rtg/variant/resources/complexhypercomplextemplate.fa"), new File(dir, "template"));
      final File sam = new File(dir, "sam.sam.gz");
      VariantTestUtils.bgzipAndIndex(FileHelper.resourceToString("com/rtg/variant/resources/complexhypercomplex.sam"), sam);
      final File output = new File(dir, "variant_out");
      final MainResult r = MainResult.run(new SingletonCli(), "-t", template.getPath(),
        "-o", output.getPath(), "--keep-duplicates", "--Xhyper-complex-length", "21",
        sam.getPath(), "--" + AbstractMultisampleCli.NO_CALIBRATION);
      assertEquals(r.err(), 0, r.rc());

      final String result = StringUtils.grep(FileHelper.gzFileToString(new File(output, "snps.vcf.gz")), "^[^#]").trim();
      //turn on following for singleton
      TestUtils.containsAll(result, "chr20 114 . C A 86.0 RX".replaceAll(" ", "\t"));
      final String bed = StringUtils.grep(FileHelper.gzFileToString(new File(output, "regions.bed.gz")), "^[^#]").trim();
      assertEquals("chr20 93 115 hyper-complex".replaceAll(" ", "\t"), bed);
    }
  }

  public void testComplexAmbiguousSuppression() throws Exception {
    try (final TestDirectory dir = new TestDirectory("test")) {
      final File template = ReaderTestUtils.getDNADir(FileHelper.resourceToString("com/rtg/variant/resources/complexambigsuppressiontemplate.fa"), new File(dir, "template"));
      final File sam = new File(dir, "sam.sam.gz");
      VariantTestUtils.bgzipAndIndex(FileHelper.resourceToString("com/rtg/variant/resources/complexambigsuppression.sam"), sam);
      final File output = new File(dir, "variant_out");
      final MainResult r = MainResult.run(new SingletonCli(), "-t", template.getPath(), "-o", output.getPath(), "--keep-duplicates", sam.getPath(), "--filter-ambiguity", "8%", "--" + AbstractMultisampleCli.NO_CALIBRATION);
      assertEquals(r.err(), 0, r.rc());

      final String result = StringUtils.grep(FileHelper.gzFileToString(new File(output, "snps.vcf.gz")), "^[^#]").trim();
      mNano.check("variantnanotest-complexambsuppression.vcf", result, false);
    }
  }

  // 12 3456 7890 1234 5678 90
  // aa atcg actg gtca gcta gg
  // atcg accg
  // cg accg gt
  // cg accg gt
  // cg accg gt
  // cg accg gt
  // accg gtca
  // acg gtcagc
  // Test handling of coverage threshold
  private static final String SAM_COMPLEX_COVERAGE_2 = "" + SAMHEADER1
      + "00a" + TAB + "0" + TAB + "g1" + TAB + "1" + TAB + "255" + TAB + "7=1X" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "AAATCGAT" + TAB + "````````" + TAB + "AS:i:1" + LS
      + SAM_COMPLEX_BODY
      + "3a" + TAB + "0" + TAB + "g1" + TAB + "7" + TAB + "255" + TAB + "1=1X1D6=" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ATGGTCAG" + TAB + "````````" + TAB + "AS:i:1" + LS
      + "3a" + TAB + "0" + TAB + "g1" + TAB + "7" + TAB + "255" + TAB + "1=1D1X6=" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ACGGTCAG" + TAB + "````````" + TAB + "AS:i:1" + LS
      + "4c" + TAB + "0" + TAB + "g1" + TAB + "9" + TAB + "255" + TAB + "1X7=" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGGTCAGC" + TAB + "````````" + TAB + "AS:i:1" + LS
      ;

  public void test56ComplexPartialCoverageThreshold2() throws Exception {
    final String[] args0 = {"--filter-depth", "8"};
    check(REF_SEQS, SAM_COMPLEX_COVERAGE_2, 56, null, 0, 88L/*regression*/, args0);
  }
  public void test56ComplexPartialCoverageThresholdUnderCoverage() throws Exception {
    final String[] args0 = {"--filter-depth", "9"};
    check(REF_SEQS, SAM_COMPLEX_COVERAGE_2, 57, null, 0, 88L/*regression*/, args0);
  }

  //site specific population priors
  public void test58() throws Exception {
    final String sam = FileHelper.resourceToString("com/rtg/variant/resources/vcftest.sam");
    final String tmplfa = FileHelper.resourceToString("com/rtg/variant/resources/vcftesttmpl.fa");
    try (final TestDirectory dir = new TestDirectory("testcheck")) {
      final File popFile = new File(dir, "pop.vcf.gz");
      BgzipFileHelper.resourceToBgzipFile("com/rtg/variant/resources/pop58.vcf", popFile);
      final File alleleCountFile = new File(dir, "allele.ac");
      final AlleleCountsFileConverter blah = new AlleleCountsFileConverter();
      blah.convert(popFile, alleleCountFile);

      check(tmplfa, sam, 58, null, "", 0, 2, null, true, 172L/*regression*/, "--population-priors", alleleCountFile.getPath());
    }
  }

  //test 3-long insert with only reads with soft clips at the start
  public void test59() throws Exception {
    final String sam = FileHelper.resourceToString("com/rtg/variant/resources/vcfcomplexsoftcliptest.sam");
    final String tmplfa = FileHelper.resourceToString("com/rtg/variant/resources/vcftesttmpl.fa");
    check(tmplfa, sam, 59, null, "", 0, 2, null, true, 120L/*regression*/);
  }
  //test a mismatch with only reads with soft clips at the start
  public void test60() throws Exception {
    final String sam = FileHelper.resourceToString("com/rtg/variant/resources/vcfsnpsoftcliptest.sam");
    final String tmplfa = FileHelper.resourceToString("com/rtg/variant/resources/vcftesttmpl.fa");
    check(tmplfa, sam, 60, null, "", 0, 2, null, true, 90L/*regression*/);
  }

  //test what happens with a really long split read generated del and some messyness from the mated file.
  public void test62() throws Exception {
    final String tmpl = ">g1\n"
        + "gtcacctttaaaagtttattgatcttttgtgacatgcacgtgggttcccagtagcaagaaactaaagggtcgcaggccggtttctgctaatttctttaa"
        + "TTCCAAGACAGTCTCAAATATTTTCTTATTAACTTCCTGGAGGGAGGCTTATCATTCTCTCTTTTGGATGATTCTAAGTACCAGCTAAAATACAGCTATCA"
        + "TTCATTTTCCTTGATTTGGG"        //deleted section
        + "AGCCTAATTTCTTTAATTTAGTATGCAAGAAAACCAATTTGGAAATATCAACTGTTTTGGAAACCTTAGACCTAGGTCA"
        + "TCCTTAGTAAGATCTTCCCATTTATATAAATACTTGCAAGTAGTAGTGCCATAATTACCAAACATAAAGCCAACTGAGATGCCCAAAGGGGGCCACTCTC"
        + "\n";

    final String sam = FileHelper.resourceToString("com/rtg/variant/resources/longindelmessy.sam");
    check(tmpl, sam, 62, "", "", 0, 2, null, true, 1900L/*regression*/);
  }

  public void testComplexSiteSpecificPriorsVcf() throws Exception {
    complexSiteSpecificPriorsCheck(true);
  }
  public void testComplexSiteSpecificPriorsAlleleCounts() throws Exception {
    complexSiteSpecificPriorsCheck(false);
  }

  public void complexSiteSpecificPriorsCheck(boolean vcf) throws Exception { //users tags 63 and 64
    //complex priors are only used to extend complex region
    final String sam = FileHelper.resourceToString("com/rtg/variant/resources/complexssp2.sam");
    final String tmplfa = FileHelper.resourceToString("com/rtg/variant/resources/complexssp.fa");
    try (final TestDirectory dir = new TestDirectory("testcheck")) {
      final File popFile = new File(dir, "pop.vcf.gz");
      BgzipFileHelper.resourceToBgzipFile("com/rtg/variant/resources/complexssp_pop.vcf", popFile);

      final File fileToUse;
      if (vcf) {
        final TabixIndexer ti = new TabixIndexer(popFile);
        ti.saveVcfIndex();
        fileToUse = popFile;
      } else {
        fileToUse = new File(dir, "allele.ac");

        final AlleleCountsFileConverter blah = new AlleleCountsFileConverter();
        blah.convert(popFile, fileToUse);
      }
      check(tmplfa, sam, 63, null, "", 0, 5, null, true, 1111L/*regression*/);
      check(tmplfa, sam, 64, null, "", 0, 5, null, true, 1111L/*regression*/, "--population-priors", fileToUse.getPath());
    }
  }

  public void test67() throws Exception {
    final String[] args0 = {"-a"};
    check(REF_SEQS67, SAM1_II, 67, null, 0, SAM1_I_LENGTH, args0);
  }


  private static final String SAM68 = ""
      + SAMHEADER11
      + "1a" + TAB + "0" + TAB + "g1" + TAB + "5" + TAB + "255" + TAB + "6M1D3M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACTGTCA" + TAB + "`````````" + TAB + "AS:i:1" + LS
      + "1b" + TAB + "0" + TAB + "g1" + TAB + "5" + TAB + "255" + TAB + "6M1D3M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACTGTCA" + TAB + "`````````" + TAB + "AS:i:1" + LS
      + "1c" + TAB + "0" + TAB + "g1" + TAB + "5" + TAB + "255" + TAB + "6M1D3M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACTGTCA" + TAB + "`````````" + TAB + "AS:i:1" + LS
      + "1d" + TAB + "0" + TAB + "g1" + TAB + "5" + TAB + "255" + TAB + "6M1D3M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACTGTCA" + TAB + "`````````" + TAB + "AS:i:1" + LS
      + "1e" + TAB + "0" + TAB + "g1" + TAB + "5" + TAB + "255" + TAB + "6M1D3M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "CGACTGTCA" + TAB + "`````````" + TAB + "AS:i:1" + LS
      + "1f" + TAB + "67" + TAB + "g1" + TAB + "5" + TAB + "255" + TAB + "6M1D3M" + TAB + "=" + TAB + "10" + TAB + "15" + TAB + "CGACTGTCA" + TAB + "`````````" + TAB + "AS:i:1" + LS
      + "1g" + TAB + "67" + TAB + "g1" + TAB + "5" + TAB + "255" + TAB + "6M1D3M" + TAB + "=" + TAB + "10" + TAB + "15" + TAB + "CGACTGACA" + TAB + "`````````" + TAB + "AS:i:2" + LS
      + "1f" + TAB + "131" + TAB + "g1" + TAB + "10" + TAB + "255" + TAB + "6M1D3M" + TAB + "=" + TAB + "5" + TAB + "-15" + TAB + "GTCAGCTAG" + TAB + "`````````" + TAB + "AS:i:1" + LS
      + "1g" + TAB + "131" + TAB + "g1" + TAB + "10" + TAB + "255" + TAB + "6M1D3M" + TAB + "=" + TAB + "5" + TAB + "-15" + TAB + "GTCAGCTAA" + TAB + "`````````" + TAB + "AS:i:2" + LS;

  //test duplicate detection and logging
  public void test68() throws Exception {
    check(REF_SEQS, SAM68, 68, null, "2 records skipped due to duplicate detection", 0, 0, null, true, 9 * 9, "-a");
  }

  //test AVR and thresholding
  public void test69() throws Exception {
    try (TestDirectory test = new TestDirectory()) {
      final File avr = new File(test, "avr");
      FileHelper.resourceToFile("com/rtg/variant/avr/resources/default.avr", avr);
      check(REF_SEQS, SAM68, 69, null, "2 records skipped due to duplicate detection", 0, 0, null, true, 9 * 9, "-a", "--avr-model", avr.getPath(), "--min-avr-score", "0.3");
    }
  }

  //test calling a bunch of Ns
  public void test70() throws Exception {
    final String tmpl = ">simulatedSequence1\n"
      + "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
      + "TGAGTGTCGTGCTTGCGGCATCCCAGTAATCATTCTGTAATATTGAAAACAACGCATAGAAAAAAAGTAGACTGTCCTGCGTGTGCTTACTAACTATCAA"
      + "CCTTGGGNGAATTAAGAGTAGCCAAGTGGCAAGACCTATGTGTACGAACCAGGGTTCATCAAGACGCGTGTATGGAGCTCTACTGTTACCCGAGCAGTTG"
      + "CGTACNNNNTGAGAGTCCCAGAGAGGTCCGTGGGAAATTTTGTGCAGTGGAGTATATGCTTTCTGAACTCCGTTAGGTCCTCTCGCATATTAACGTAGGG"
      + "CTAATGATAGTGCANCACAGTGCTACCTATAATTCTCACGAGCTGCACCCCCCGAAAGGATGCCACTCATGCCGTCCCTGGCGGGGAGCAGTAAAATAGT"
      + "TACATCACCGCTACCCAAGTGTCAAACGTCTTTGGCCGGGAGCCAGGCCCACGTGATGGACCTTCATGCNNNGTTTCCGTCGAGGCCGAGTCCCTATTTA"
      + "CTATCGGGGGAATTCATACCCTGAGATAAAACACCTTCCTGATACGGGCGTCGACCGCAATCCCTCCCCGAGTAATCCTTGGTTGAGAATCACCCCGACT"
      + "NNNACAGAGCTCGTCCGGATTCAAAGTAGGATTCCCGCATCGGTCTGAGGCGGTCGAGCCACTCATCATGTGGAAGTGTATGAAGGACTTAAAAACTGAG"
      + "GCTACCGCATATTAGAGGGCACAACCCCGCAGCCGNNCTAATAGCAATCGAAGTCTAACCGGTAGTACTCAGGCTGAGCTCCCCCTCGCGCTTGTATGTT"
      + "TAAAGCCTTGTGGTATAGTTGCGGTTAAGCTTAGAAACAGACCCCCCCCAGTAATGTCGATCGTAGCAACTAAAACGTGTCTGGTGAGTAAAGCAGAGTC\n";

    final String sam = FileHelper.resourceToString("com/rtg/variant/resources/n-mappings.sam");
    check(tmpl, sam, 70, "", "", 0, 2, null, true, -1L);
  }

  public void test72va() throws Exception {
    final String tmpl = ">simulatedSequence1\n"
      + "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
      + "TGAGTGTCGTGCTTGCGGCATCCCAGTAATCATTCTGTAATATTGAAAACAACGCATAGAAAAAAAGTAGACTGTCCTGCGTGTGCTTACTAACTATCAA"
      + "CCTTGGGNGAATTAAGAGTAGCCAAGTGGCAAGACCTATGTGTACGAACCAGGGTTCATCAAGACGCGTGTATGGAGCTCTACTGTTACCCGAGCAGTTG"
      + "CGTACNNNNTGAGAGTCCCAGAGAGGTCCGTGGGAAATTTTGTGCAGTGGAGTATATGCTTTCTGAACTCCGTTAGGTCCTCTCGCATATTAACGTAGGG"
      + "CTAATGATAGTGCANCACAGTGCTACCTATAATTCTCACGAGCTGCACCCCCCGAAAGGATGCCACTCATGCCGTCCCTGGCGGGGAGCAGTAAAATAGT"
      + "TACATCACCGCTACCCAAGTGTCAAACGTCTTTGGCCGGGAGCCAGGCCCACGTGATGGACCTTCATGCNNNGTTTCCGTCGAGGCCGAGTCCCTATTTA"
      + "CTATCGGGGGAATTCATACCCTGAGATAAAACACCTTCCTGATACGGGCGTCGACCGCAATCCCTCCCCGAGTAATCCTTGGTTGAGAATCACCCCGACT"
      + "NNNACAGAGCTCGTCCGGATTCAAAGTAGGATTCCCGCATCGGTCTGAGGCGGTCGAGCCACTCATCATGTGGAAGTGTATGAAGGACTTAAAAACTGAG"
      + "GCTACCGCATATTAGAGGGCACAACCCCGCAGCCGNNCTAATAGCAATCGAAGTCTAACCGGTAGTACTCAGGCTGAGCTCCCCCTCGCGCTTGTATGTT"
      + "TAAAGCCTTGTGGTATAGTTGCGGTTAAGCTTAGAAACAGACCCCCCCCAGTAATGTCGATCGTAGCAACTAAAACGTGTCTGGTGAGTAAAGCAGAGTC\n";

    final String sam = FileHelper.resourceToString("com/rtg/variant/resources/n-mappings.sam");
    check(tmpl, sam, 72, "", "", 0, 2, null, true, -1L, "--min-variant-allelic-depth", "1", "--min-variant-allelic-fraction", "0.03");
  }


  // Test situation where N's on ref and evicence is shorter sequence of AAAAA.  See Bug#1591
  public void testIndelWithUnknownRef() throws Exception {
    final String tmpl = ">ref\n"
    + "AGTTAAAGAGTGAAACCCTGATAGTCTTACCCCAAGGCCAAAGTCCTATTTTATTATTTTTATATTCTTACTATATATTATACAAATCTTCATTGCAA"
    + "GTTTNNNNNNNNNNNNNNNNNNNNAAGTAAAAACATAAGAAATCTAATTTTTGTATATAAAAGCTGTAAACTAAATTATATATATACACATACATACA"
    + "TACGTGTGTGTGTGTATATATATATACATATATAACCTATGGATTAGGAAAATTTATTGCTTCAACAAACTAAGGGGATTACTTCCCATAAAATTAGT\n";
    final String sam = FileHelper.resourceToString("com/rtg/variant/resources/n2-mappings.sam");
    check(tmpl, sam, 71, "", "", 0, 2, null, true, -1L);
  }

  //See bugzilla #1524, this test replicates the problem that shows up in 2.7 but has been fixed at some point and so now works
  public void testBug1524() throws Exception {
    try (final TestDirectory dir = new TestDirectory()) {
      final File sam = FileHelper.resourceToFile("com/rtg/variant/resources/bug1524.sam.gz", new File(dir, "bug1524.sam.gz"));
      final File genome = ReaderTestUtils.getDNADir(FileHelper.resourceToString("com/rtg/variant/resources/bug1524.fasta"), new File(dir, "genome"));
      new TabixIndexer(sam).saveSamIndex();
      final File output = new File(dir, "variant_out");
      final MainResult r = MainResult.run(new SingletonCli(), "-Z", "-t", genome.getPath(), "-o", output.getPath(), sam.getPath(), "--" + AbstractMultisampleCli.NO_CALIBRATION);
      assertEquals(r.err(), 0, r.rc());
      TestUtils.containsAll(r.err(), "SAM record is invalid", "1 records skipped because of SAM format problems");
      mNano.check("bug1524.txt", r.out());
      final String result = FileUtils.fileToString(new File(output, VariantParams.VCF_OUT_SUFFIX));
      final String actualFixed = TestUtils.sanitizeVcfHeader(result);
      mNano.check("bug1524.vcf", actualFixed, false);
    }
  }

}
