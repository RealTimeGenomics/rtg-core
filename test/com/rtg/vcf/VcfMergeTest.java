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
package com.rtg.vcf;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

import com.rtg.tabix.TabixIndexer;
import com.rtg.tabix.UnindexableDataException;
import com.rtg.util.intervals.RegionRestriction;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.BgzipFileHelper;
import com.rtg.util.test.FileHelper;
import com.rtg.util.test.NanoRegression;
import com.rtg.variant.VariantStatistics;
import com.rtg.vcf.header.MetaType;
import com.rtg.vcf.header.VcfHeader;
import com.rtg.vcf.header.VcfNumber;

import junit.framework.TestCase;

/**
 * Test class
 */
public class VcfMergeTest extends TestCase {

  private NanoRegression mNano;

  @Override
  public void setUp() {
    Diagnostic.setLogStream();
    mNano = new NanoRegression(VcfMergeTest.class);
  }

  @Override
  public void tearDown() throws Exception {
    try {
      mNano.finish();
    } finally {
      mNano = null;
    }
  }

  private static final VcfHeader MINIMAL_HEADER;
  private static final VcfHeader HEADER2;
  static {
    MINIMAL_HEADER = new VcfHeader();
    MINIMAL_HEADER.setVersionValue("VCFv4.1");
    MINIMAL_HEADER.addFormatField("GT", MetaType.STRING, new VcfNumber("1"), "Genotype");
    MINIMAL_HEADER.addLine(VcfHeader.META_STRING + "contig=<ID=\"chr2\",length=10000>");
    MINIMAL_HEADER.addLine(VcfHeader.META_STRING + "contig=<ID=\"chr3\",length=10000>");
    MINIMAL_HEADER.addLine(VcfHeader.META_STRING + "contig=<ID=\"1\",length=10000000>");
    MINIMAL_HEADER.addSampleName("%1$s");

    HEADER2 = new VcfHeader();
    HEADER2.setVersionValue("VCFv4.1");
    HEADER2.addFormatField("GT", MetaType.STRING, new VcfNumber("1"), "Genotype");
    HEADER2.addLine(VcfHeader.META_STRING + "INFO=<ID=XRX,Number=0,Type=Flag,Description=\"RTG variant was called using complex caller\">");
    HEADER2.addLine(VcfHeader.META_STRING + "INFO=<ID=RCE,Number=0,Type=Flag,Description=\"RTG variant is equivalent to the previous variant\">");
    HEADER2.addLine(VcfHeader.META_STRING + "INFO=<ID=CT,Number=1,Type=Integer,Description=\"Coverage threshold that was applied\">");
    HEADER2.addLine(VcfHeader.META_STRING + "contig=<ID=\"chr2\",length=10000>");
    HEADER2.addLine(VcfHeader.META_STRING + "contig=<ID=\"chr3\",length=10000>");
    HEADER2.addSampleName("%1$s");
  }

  private static final String VCF_RECORD_FORMAT = "%1$s\t%2$d\t.\t%3$s\t%4$s\t.\tPASS\t%7$s\tGT\t%5$d/%6$d";

  private static String createRecord(String chrom, int pos, String ref, int gt1, int gt2, String info, String... alts) {
    final String altsStr = StringUtils.implode(alts, ",", false);
    return String.format(VCF_RECORD_FORMAT, chrom, pos, ref, altsStr, gt1, gt2, info);
  }
  private static final String FILE1 = String.format(MINIMAL_HEADER.toString(), "sample1")
            + createRecord("chr2", 100, "A", 1, 1, ".", "C") + "\n"
            + createRecord("chr2", 200, "A", 1, 1, ".", "C") + "\n"
            + createRecord("chr3", 300, "A", 1, 1, ".", "C") + "\n"
            + createRecord("chr3", 500, "A", 1, 1, ".", "C") + "\n";
  private static final String FILE2 = String.format(MINIMAL_HEADER.toString(), "sample2")
            + createRecord("chr3", 200, "A", 1, 1, ".", "G") + "\n"
            + createRecord("chr3", 500, "A", 1, 1, ".", "G") + "\n"
            + createRecord("chr2", 100, "A", 1, 1, ".", "G") + "\n"
            + createRecord("chr2", 900, "A", 1, 1, ".", "G") + "\n";
  private static final String FILE3 = String.format(HEADER2.toString(), "sample3")
      + createRecord("chr3", 200, "A", 1, 1, ".", "G") + "\n"
      + createRecord("chr3", 500, "A", 1, 1, "XRX", "G") + "\n"
      + createRecord("chr2", 100, "A", 1, 1, "RCE", "G") + "\n"
      + createRecord("chr2", 900, "A", 1, 1, "CT=123", "G") + "\n";
  public void testZipper() throws Exception {
    try (final TestDirectory dir = new TestDirectory("vcfmerge")) {
      final File f1 = BgzipFileHelper.bytesToBgzipFile(FILE1.getBytes(), new File(dir, "file1.vcf.gz"));
      final File f2 = BgzipFileHelper.bytesToBgzipFile(FILE2.getBytes(), new File(dir, "file2.vcf.gz"));
      new TabixIndexer(f1, TabixIndexer.indexFileName(f1)).saveVcfIndex();
      new TabixIndexer(f2, TabixIndexer.indexFileName(f2)).saveVcfIndex();

      final int[] positions = {100, 200, 900, 200, 300, 500};
      final String[] templates = {"chr2", "chr2", "chr2", "chr3", "chr3", "chr3"};
      final String[][] alts = {new String[] {"C", "G"}, new String[] {"C"}, new String[] {"G"},
        new String[] {"G"}, new String[] {"C"}, new String[] {"C", "G"}};
      final String[][] samples = {new String[] {"sample1", "sample2"}, new String[] {"sample1"}, new String[] {"sample2"},
        new String[] {"sample2"}, new String[] {"sample1"}, new String[] {"sample1", "sample2"}};

      VcfMerge.VcfPositionZipper zipper = new VcfMerge.VcfPositionZipper(null, f1, f2);
      try {
        final VcfMerge.ZipperCallback callback = new ZipperCallbackTester(positions, templates, alts, samples);
        int count = 0;
        while (zipper.hasNextPosition()) {
          zipper.nextPosition(callback);
          count++;
        }
        assertEquals(positions.length, count);
      } finally {
        zipper.close();
      }
      zipper = new VcfMerge.VcfPositionZipper(new RegionRestriction("chr3:300-500"), f1, f2);
      try {
        final VcfMerge.ZipperCallback callback = new ZipperCallbackTester(
                Arrays.copyOfRange(positions, 4, 6),
                Arrays.copyOfRange(templates, 4, 6),
                Arrays.copyOfRange(alts, 4, 6),
                Arrays.copyOfRange(samples, 4, 6));
        int count = 0;
        while (zipper.hasNextPosition()) {
          zipper.nextPosition(callback);
          count++;
        }
        assertEquals(2, count);
      } finally {
        zipper.close();
      }
    }
  }

  public void testMerge() throws Exception {
    try (final TestDirectory dir = new TestDirectory("vcfmerge")) {
      final File f1 = BgzipFileHelper.bytesToBgzipFile(FILE1.getBytes(), new File(dir, "file1.vcf.gz"));
      final File f2 = BgzipFileHelper.bytesToBgzipFile(FILE2.getBytes(), new File(dir, "file2.vcf.gz"));
      new TabixIndexer(f1, TabixIndexer.indexFileName(f1)).saveVcfIndex();
      new TabixIndexer(f2, TabixIndexer.indexFileName(f2)).saveVcfIndex();
      final File output = new File(dir, "out.vcf");
      VcfMerge.mergeVcfFiles(null, output, false, false, null, new String[]{"##extraline=foo", "##extraline2=bar" }, null, false, f1, f2);
      String actual = FileUtils.fileToString(output);
      actual = StringUtils.grepMinusV(actual, "^##(RUN-ID)|(CL)").replaceAll("[\r\n]+", "\n");
      mNano.check("vcfmerge_testMerge.vcf", actual, false);
    }
  }

  public void testMergeFileOrderReversed() throws Exception {
    try (final TestDirectory dir = new TestDirectory("vcfmerge")) {
      final File f2 = BgzipFileHelper.bytesToBgzipFile(FILE1.getBytes(), new File(dir, "file1.vcf.gz"));
      final File f1 = BgzipFileHelper.bytesToBgzipFile(FILE2.getBytes(), new File(dir, "file2.vcf.gz"));
      new TabixIndexer(f1, TabixIndexer.indexFileName(f1)).saveVcfIndex();
      new TabixIndexer(f2, TabixIndexer.indexFileName(f2)).saveVcfIndex();
      final File output = new File(dir, "out.vcf");
      final VariantStatistics stats = new VariantStatistics(null);
      VcfMerge.mergeVcfFiles(null, output, true, true, stats, new String[]{"##extraline=foo", "##extraline2=bar" }, null, false, f1, f2);
      String actual = FileHelper.gzFileToString(new File(dir, "out.vcf.gz"));
      actual = StringUtils.grepMinusV(actual, "^##(RUN-ID)|(CL)").replaceAll("[\r\n]+", "\n");
      assertTrue(new File(dir, "out.vcf.gz.tbi").isFile());
      mNano.check("vcfmerge_testMerge_1.vcf", actual, false);
      final MemoryPrintStream ps = new MemoryPrintStream();
      stats.printStatistics(ps.outputStream());
      mNano.check("vcfmerge_stats.txt", ps.toString());
    }
  }


  public void testSingleFileMerge() throws Exception {
    try (final TestDirectory dir = new TestDirectory("vcfmerge")) {
      final File f1 = BgzipFileHelper.bytesToBgzipFile(FILE3.getBytes(), new File(dir, "file1.vcf.gz"));
      new TabixIndexer(f1, TabixIndexer.indexFileName(f1)).saveVcfIndex();
      final File output = new File(dir, "out.vcf");
      VcfMerge.mergeVcfFiles(null, output, false, false, null, new String[]{"##extraline=foo", "##extraline2=bar" }, null, false, f1);
      String actual = FileUtils.fileToString(output);
      actual = StringUtils.grepMinusV(actual, "^##(RUN-ID)|(CL)").replaceAll("[\r\n]+", "\n");
      mNano.check("vcfmerge_testSingle.vcf", actual, false);
    }
  }


  private static class ZipperCallbackTester implements VcfMerge.ZipperCallback {
    private int mI = 0;
    private final int[] mPositions;
    private final String[] mTemplates;
    private final String[][] mAlts;
    private final String[][] mSamples;

    public ZipperCallbackTester(int[] positions, String[] templates, String[][] alts, String[][] samples) {
      mPositions = positions;
      mTemplates = templates;
      mAlts = alts;
      mSamples = samples;
    }

    @Override
    public void vcfAtPosition(VcfRecord[] records, VcfHeader[] headers) {
      for (int i = 0; i < records.length; i++) {
        assertEquals("" + mI, mPositions[mI], records[i].getOneBasedStart());
        assertEquals(mTemplates[mI], records[i].getSequenceName());
        assertEquals(mAlts[mI][i], records[i].getAltCalls().get(0));
        assertTrue(headers[i].toString().contains(mSamples[mI][i]));
      }
      mI++;
    }
  }

  public void testSamePositionDifferentRef() throws Exception {
    try (final TestDirectory dir = new TestDirectory("vcfmerge")) {
      final String string1 = String.format(MINIMAL_HEADER.toString(), "NA12878")
          + "1       20200602        .       C       T       575.1   PASS    XRX     GT:DP:RE:GQ     0/1:42:4.891:56\n".replaceAll(" +", "\t")
          + "1       20200605        .       T       C       575.1   PASS    XRX     GT:DP:RE:GQ     0/1:42:4.891:56\n".replaceAll(" +", "\t")
          + "1       20200635        .       C       T       575.1   PASS    XRX     GT:DP:RE:GQ     1/1:42:4.891:56\n".replaceAll(" +", "\t");

      final String string2 = String.format(MINIMAL_HEADER.toString(), "NA12891")
          + "1 20200613        .       ACCACCATCAC     ACCACCAT,ATCACCATCAC    895.7   PASS    XRX     GT:DP:RE:GQ     1/2:54:2.680:43\n".replaceAll(" +", "\t")
          + "1       20200635        .       C       CCAT    601.1   PASS    XRX     GT:DP:RE:GQ     0/1:46:0.028:601\n".replaceAll(" +", "\t");

      final String string3 = String.format(MINIMAL_HEADER.toString(), "NA12892")
          + "1       20200602        .       CCATCACCACCACCACCATCACCACCACCACCAC      CCATCACCACCACCACCAT,CTATCACCACCATCACCATCACCACCACCACCACCAT       2562.9  PASS    XRX     GT:DP:RE:GQ     1/2:57:2.515:371\n".replaceAll(" +", "\t");

      final File f1 = BgzipFileHelper.bytesToBgzipFile(string1.getBytes(), new File(dir, "file1.vcf.gz"));
      new TabixIndexer(f1, TabixIndexer.indexFileName(f1)).saveVcfIndex();
      final File f2 = BgzipFileHelper.bytesToBgzipFile(string2.getBytes(), new File(dir, "file2.vcf.gz"));
      new TabixIndexer(f2, TabixIndexer.indexFileName(f2)).saveVcfIndex();
      final File f3 = BgzipFileHelper.bytesToBgzipFile(string3.getBytes(), new File(dir, "file3.vcf.gz"));
      new TabixIndexer(f3, TabixIndexer.indexFileName(f3)).saveVcfIndex();
      final File output = new File(dir, "out.vcf");
      VcfMerge.mergeVcfFiles(null, output, false, false, null, new String[]{}, null, false, f1, f2, f3);
      String actual = FileUtils.fileToString(output);
      actual = StringUtils.grepMinusV(actual, "^##(RUN-ID)|(CL)").replaceAll("[\r\n]+", "\n");
      mNano.check("vcfmerge_testSamePosDiffRef.vcf", actual, false);
    }
  }

  public void testSamePosition() throws Exception {
    try (final TestDirectory dir = new TestDirectory("vcfmerge")) {
      final String string1 = String.format(MINIMAL_HEADER.toString(), "NA12878")
              + "1       20200602        .       C       T       575.1   PASS    XRX     GT:DP:RE:GQ     0/1:42:4.891:56\n".replaceAll(" +", "\t");
      final String string2 = String.format(MINIMAL_HEADER.toString(), "NA12891")
              + "1       20200602        .       CA       C       575.1   PASS    XRX     GT:DP:RE:GQ     0/1:42:4.891:56\n".replaceAll(" +", "\t");
      final String string3 = String.format(MINIMAL_HEADER.toString(), "NA12892")
              + "1       20200602        .       C       G       575.1   PASS    XRX     GT:DP:RE:GQ     0/1:42:4.891:56\n".replaceAll(" +", "\t");
      final File f1 = BgzipFileHelper.bytesToBgzipFile(string1.getBytes(), new File(dir, "file1.vcf.gz"));
      new TabixIndexer(f1, TabixIndexer.indexFileName(f1)).saveVcfIndex();
      final File f2 = BgzipFileHelper.bytesToBgzipFile(string2.getBytes(), new File(dir, "file2.vcf.gz"));
      new TabixIndexer(f2, TabixIndexer.indexFileName(f2)).saveVcfIndex();
      final File f3 = BgzipFileHelper.bytesToBgzipFile(string3.getBytes(), new File(dir, "file3.vcf.gz"));
      new TabixIndexer(f3, TabixIndexer.indexFileName(f3)).saveVcfIndex();
      final File output = new File(dir, "out.vcf.gz");
      VcfMerge.mergeVcfFiles(null, output, true, true, null, new String[]{}, null, false, f2, f3);
      final File output2 = new File(dir, "out2.vcf");
      VcfMerge.mergeVcfFiles(null, output2, false, false, null, new String[]{}, null, false, f1, output);
      String actual = FileUtils.fileToString(output2);
      actual = StringUtils.grepMinusV(actual, "^##(RUN-ID)|(CL)").replaceAll("[\r\n]+", "\n");
      mNano.check("vcfmerge_testSamePos.vcf", actual, false);
    }

  }

  public void testMergingVariableNumberFormats() throws IOException, UnindexableDataException {
    try (final TestDirectory dir = new TestDirectory()) {
      final File snpsA = BgzipFileHelper.bytesToBgzipFile(FileHelper.resourceToString("com/rtg/vcf/resources/snpsA.vcf").getBytes(), new File(dir, "snpsA.vcf.gz"));
      new TabixIndexer(snpsA, TabixIndexer.indexFileName(snpsA)).saveVcfIndex();
      final File snpsB = BgzipFileHelper.bytesToBgzipFile(FileHelper.resourceToString("com/rtg/vcf/resources/snpsB.vcf").getBytes(), new File(dir, "snpsB.vcf.gz"));
      new TabixIndexer(snpsB, TabixIndexer.indexFileName(snpsB)).saveVcfIndex();
      final File outFile = new File(dir, "test.vcf");
      VcfMerge.mergeVcfFiles(null, outFile, false, false, null, new String[]{}, null, true, snpsA, snpsB);
      String actual = FileUtils.fileToString(outFile);
      actual = StringUtils.grepMinusV(actual, "^##(RUN-ID)|(CL)").replaceAll("[\r\n]+", "\n");
      mNano.check("mergedAB.vcf", actual, false);
    }
  }

  public void testMergeDrop() throws Exception {
    try (final TestDirectory dir = new TestDirectory("vcfmerge")) {
      final File snpsA = BgzipFileHelper.bytesToBgzipFile(FileHelper.resourceToString("com/rtg/vcf/resources/vcfmerge-na12889.vcf").getBytes(), new File(dir, "na12889.vcf.gz"));
      new TabixIndexer(snpsA, TabixIndexer.indexFileName(snpsA)).saveVcfIndex();
      final File snpsB = BgzipFileHelper.bytesToBgzipFile(FileHelper.resourceToString("com/rtg/vcf/resources/vcfmerge-na12880.vcf").getBytes(), new File(dir, "na12880.vcf.gz"));
      new TabixIndexer(snpsB, TabixIndexer.indexFileName(snpsB)).saveVcfIndex();
      final File output = new File(dir, "out.vcf");
      final VariantStatistics stats = new VariantStatistics(null);
      VcfMerge.mergeVcfFiles(null, output, true, true, stats, new String[]{"##extraline=foo", "##extraline2=bar" }, null, false, snpsA, snpsB);
      String actual = FileHelper.gzFileToString(new File(dir, "out.vcf.gz"));
      actual = StringUtils.grepMinusV(actual, "^##(RUN-ID)|(CL)").replaceAll("[\r\n]+", "\n");
      assertTrue(new File(dir, "out.vcf.gz.tbi").isFile());
      mNano.check("vcfmerge_testMerge_drop.vcf", actual, false);
      final MemoryPrintStream ps = new MemoryPrintStream();
      stats.printStatistics(ps.outputStream());
      mNano.check("vcfmerge_stats_drop.txt", ps.toString());
    }
  }

  public void testMergePreserve() throws Exception {
    try (final TestDirectory dir = new TestDirectory("vcfmerge")) {
      final File snpsA = BgzipFileHelper.bytesToBgzipFile(FileHelper.resourceToString("com/rtg/vcf/resources/vcfmerge-na12889.vcf").getBytes(), new File(dir, "na12889.vcf.gz"));
      new TabixIndexer(snpsA, TabixIndexer.indexFileName(snpsA)).saveVcfIndex();
      final File snpsB = BgzipFileHelper.bytesToBgzipFile(FileHelper.resourceToString("com/rtg/vcf/resources/vcfmerge-na12880.vcf").getBytes(), new File(dir, "na12880.vcf.gz"));
      new TabixIndexer(snpsB, TabixIndexer.indexFileName(snpsB)).saveVcfIndex();
      final File output = new File(dir, "out.vcf");
      final VariantStatistics stats = new VariantStatistics(null);
      VcfMerge.mergeVcfFiles(null, output, true, true, stats, new String[]{"##extraline=foo", "##extraline2=bar" }, null, true, snpsA, snpsB);
      String actual = FileHelper.gzFileToString(new File(dir, "out.vcf.gz"));
      actual = StringUtils.grepMinusV(actual, "^##(RUN-ID)|(CL)").replaceAll("[\r\n]+", "\n");
      assertTrue(new File(dir, "out.vcf.gz.tbi").isFile());
      mNano.check("vcfmerge_testMerge_preserve.vcf", actual, false);
      final MemoryPrintStream ps = new MemoryPrintStream();
      stats.printStatistics(ps.outputStream());
      mNano.check("vcfmerge_stats_preserve.txt", ps.toString());
    }
  }

}
