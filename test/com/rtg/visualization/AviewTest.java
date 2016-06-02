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

package com.rtg.visualization;

import static com.rtg.util.StringUtils.LS;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.launcher.globals.GlobalFlags;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SdfId;
import com.rtg.sam.Sam2Bam;
import com.rtg.util.TestUtils;
import com.rtg.util.Utils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;
import com.rtg.vcf.header.VcfHeader;

/**
 */
public class AviewTest extends AbstractCliTest {

  private static final String T = "\t";
  static final String ALIGNMENTS_BAM_FILE_NAME = "alignments.bam";
  static final String ALIGNMENTS_SAM_FILE_NAME = "alignments.sam";
  static final String TEMPLATE_FILE_NAME = "template";
  static final String SNPS_FILE_NAME = "snps.vcf";
  static final String BED_FILE_NAME = "snps.bed";

  /** Reference Sequence **/
  static final String REF_SEQS = "" + ">g1" + LS + "aa" + "atcg" + "actg" + "gtca" + "gcta" + "gg" + LS;

  static final String SAM = "" + "@HD" + T + "VN:1.0" + T + "SO:coordinate" + LS
  + "@SQ" + T + "SN:g1" + T + "LN:20" + LS
  + "0" + T + "0" + T + "g1" + T + "3" + T + "255" + T + "8M" + T + "*" + T + "0" + T + "0" + T + "ATCGACTG" + T + "&'(`````" + T + "AS:i:0" + T + "IH:i:1" + LS
  + "1" + T + "0" + T + "g1" + T + "3" + T + "255" + T + "8M" + T + "*" + T + "0" + T + "0" + T + "ATCGACTG" + T + "````````" + T + "AS:i:0" + T + "IH:i:1" + LS
  + "2" + T + "0" + T + "g1" + T + "5" + T + "255" + T + "8M" + T + "*" + T + "0" + T + "0" + T + "CGACTGGT" + T + "````````" + T + "AS:i:1" + T + "IH:i:1" + LS
  + "3" + T + "0" + T + "g1" + T + "6" + T + "255" + T + "8M" + T + "*" + T + "0" + T + "0" + T + "GACTGGTC" + T + "````````" + T + "AS:i:1" + T + "IH:i:1" + LS
  + "4" + T + "0" + T + "g1" + T + "6" + T + "255" + T + "8M" + T + "*" + T + "0" + T + "0" + T + "GACTGGTC" + T + "````````" + T + "AS:i:1" + T + "IH:i:1" + LS
  + "5" + T + "0" + T + "g1" + T + "11" + T + "255" + T + "8M" + T + "*" + T + "0" + T + "0" + T + "GTCAGCTA" + T + "````````" + T + "AS:i:1" + T + "IH:i:1" + LS;

  static final String SNP = ""
//    + "#name\tposition\ttype\treference\tprediction\tposterior\tcoverage\tcorrection\tsupport_statistics" + StringUtils.LS
      + VcfHeader.MINIMAL_HEADER + "\tSAMPLE" + LS
//    + "g1\t11\te\tG\tC:T\t9.3\t4\t0.040\tC\t2\t0.020\tT\t2\t0.020" + StringUtils.LS;
    + "g1\t11\t.\tG\tC,T\t9.3\tPASS\t.\tGT\t1/2" + LS;

  void run(String testname, final String sam, final String... extraargs) throws IOException {
    runWithErr(testname, sam, "", extraargs);
  }

  void runWithErr(String testname, final String sam, final String expectErr, final String... extraargs) throws IOException {
    Diagnostic.setLogStream();
    final File f = FileUtils.createTempDir("aview", "filetest");
    try {
      prepareData(f, sam);
      final String[] args = {
          new File(f, ALIGNMENTS_BAM_FILE_NAME).getPath(),
          "-c", new File(f, SNPS_FILE_NAME).getPath(),
          "-b", new File(f, SNPS_FILE_NAME).getPath(),
          "-t", new File(f, TEMPLATE_FILE_NAME).getPath(),
      };
      final Aview aview = new Aview();
      assertEquals("aview", aview.moduleName());
      final ByteArrayOutputStream bos = new ByteArrayOutputStream();
      final MemoryPrintStream err = new MemoryPrintStream();
      final int code2 = aview.mainInit(Utils.append(args, extraargs), bos, err.printStream());
      GlobalFlags.resetAccessedStatus();
      assertEquals(err.toString(), 0, code2);
      assertEquals(expectErr, err.toString());
      bos.flush();
      final String res = bos.toString();
      //System.err.println(res);
      mNano.check(testname, res);
      bos.close();
      err.close();

    } finally {
      deleteBrokenBam(f);
    }
  }

  static void prepareData(final File f, final String sam) throws IOException {
    final File template = new File(f, TEMPLATE_FILE_NAME);
    ReaderTestUtils.getReaderDNA(REF_SEQS, template, new SdfId(2L)).close();
    final File alignments = new File(f, ALIGNMENTS_SAM_FILE_NAME);
    FileUtils.stringToFile(sam, alignments);
    final File snps = new File(f, SNPS_FILE_NAME);
    FileUtils.stringToFile(SNP, snps);
    final File bam = new File(f, ALIGNMENTS_BAM_FILE_NAME);
    final MemoryPrintStream errStream = new MemoryPrintStream();
    final Sam2Bam bc = new Sam2Bam();
    GlobalFlags.resetAccessedStatus();
    final int code = bc.mainInit(new String[] {"-o", bam.getPath(), alignments.getPath()}, TestUtils.getNullOutputStream(), errStream.printStream());
    GlobalFlags.resetAccessedStatus();
    assertEquals(errStream.toString(), 0, code);
  }

  static void deleteBrokenBam(final File f) {
    //this ugliness is due to java NIO MappedByteByteBuffer can not release the index file immediately
    //following assertTrue can still fail on windows some times, thats why need to run System.gc
    //for more info look at http://bugs.sun.com/bugdatabase/view_bug.do?bug_id=4724038
    System.gc();

    System.runFinalization();
    assertTrue(FileHelper.deleteAll(f));
  }

  public void testNoColor() throws IOException, InterruptedException {
    run("aview-nocolor.txt", SAM, "--no-color", "--region", "g1:10+2", "--no-dots");
  }

  public void testColor() throws IOException, InterruptedException {
    run("aview-color.txt", SAM, "--region", "g1:11+4", "--no-dots", "--no-base-color");
  }

  public void testBaseColor() throws IOException, InterruptedException {
    run("aview-basecolor.txt", SAM, "--region", "g1:11+4", "--no-dots");
  }

  @Override
  protected AbstractCli getCli() {
    return new Aview();
  }

  public final void testInitFlags() {
    checkHelp("alignment SAM/BAM files",
        "reference SDF",
        "print reference line every N lines",
        "called variants",
        "baseline variants",
        "read SDF",
        "display nucleotide instead of dots",
        "print alignment cigars",
        "do not use colors",
        "print read name"
    );
  }

  static final String SAMINSERT = "" + "@HD" + T + "VN:1.0" + T + "SO:coordinate" + LS
  + "@SQ" + T + "SN:g1" + T + "LN:20" + LS
  + "0" + T + "0" + T + "g1" + T + "3" + T + "255" + T + "8=" + T + "*" + T + "0" + T + "0" + T + "ATCGACTG" + T + "&'(`````" + T + "AS:i:0" + T + "IH:i:1" + LS
  + "1" + T + "0" + T + "g1" + T + "3" + T + "255" + T + "8=" + T + "*" + T + "0" + T + "0" + T + "ATCGACTG" + T + "````````" + T + "AS:i:0" + T + "IH:i:1" + LS
  + "2" + T + "0" + T + "g1" + T + "5" + T + "255" + T + "1=1I7=" + T + "*" + T + "0" + T + "0" + T + "CAGACTGGT" + T + "`````````" + T + "AS:i:1" + T + "IH:i:1" + LS
  + "3" + T + "0" + T + "g1" + T + "6" + T + "255" + T + "8=" + T + "*" + T + "0" + T + "0" + T + "GACTGGTC" + T + "````````" + T + "AS:i:1" + T + "IH:i:1" + LS
  + "4" + T + "0" + T + "g1" + T + "6" + T + "255" + T + "8=" + T + "*" + T + "0" + T + "0" + T + "GACTGGTC" + T + "````````" + T + "AS:i:1" + T + "IH:i:1" + LS
  + "5" + T + "0" + T + "g1" + T + "11" + T + "255" + T + "8=" + T + "*" + T + "0" + T + "0" + T + "GTCAGCTA" + T + "````````" + T + "AS:i:1" + T + "IH:i:1" + LS;


  public void testInserts() throws Exception {
    run("aview-inserts.txt", SAMINSERT, "--no-color", "--region", "g1:9+3");
  }
  static final String SAMDELETE = ""
          + "@HD" + T + "VN:1.0" + T + "SO:coordinate" + LS
          + "@SQ" + T + "SN:g1" + T + "LN:20" + LS
          + "0" + T + "0" + T + "g1" + T + "8" + T + "255" + T + "5=" + T + "*" + T + "0" + T + "0" + T + "CTGGT" + T + "&'(``" + T + "AS:i:0" + T + "IH:i:1" + LS;
   static final String SNPDELETE = ""
//    + "#name\tposition\ttype\treference\tprediction\tposterior\tcoverage\tcorrection\tsupport_statistics" + StringUtils.LS
      + VcfHeader.MINIMAL_HEADER + "\tSAMPLE" + LS
//    + "g1\t11\te\tG\tC:T\t9.3\t4\t0.040\tC\t2\t0.020\tT\t2\t0.020" + StringUtils.LS;
       + "g1\t11\t.\tGTCA\tG\t9.3\tPASS\t.\tGT\t1/1" + LS;
  static final String BEDDELETE = ""
//    + chr start end label score
      + "g1\t10\t13\tcomplex-called" + LS;

  public void testDeletes() throws Exception {
    try (TestDirectory dir = new TestDirectory()) {
      final File template = new File(dir, TEMPLATE_FILE_NAME);
      ReaderTestUtils.getReaderDNA(REF_SEQS, template, new SdfId(2L)).close();

      final File alignments = new File(dir, ALIGNMENTS_SAM_FILE_NAME);
      FileUtils.stringToFile(SAMDELETE, alignments);
      final File snps = new File(dir, SNPS_FILE_NAME);
      FileUtils.stringToFile(SNPDELETE, snps);
      final File bed = new File(dir, BED_FILE_NAME);
      FileUtils.stringToFile(BEDDELETE, bed);
      final File bam = new File(dir, ALIGNMENTS_BAM_FILE_NAME);
      final MemoryPrintStream errStream = new MemoryPrintStream();
      final Sam2Bam bc = new Sam2Bam();
      final int code = bc.mainInit(new String[] {"-o", bam.getPath(), alignments.getPath()}, TestUtils.getNullOutputStream(), errStream.printStream());
      GlobalFlags.resetAccessedStatus();
      assertEquals(errStream.toString(), 0, code);
      final MemoryPrintStream result = new MemoryPrintStream();
      final int code2 = new Aview().mainInit(new String[] {"-t", template.getPath(), "--region", "g1:9", bam.getPath(), "--bed", bed.getPath() + "=mybed", "-c", snps.getPath() + "=mysnps",
"--no-color"}, result.outputStream(), errStream.printStream());
      GlobalFlags.resetAccessedStatus();
      assertEquals(errStream.toString(), 0, code2);
      mNano.check("aview-delete-exp.txt", result.toString());
      //System.err.println(errStream.toString());
    }
    run("aview-deletes.txt", SAMINSERT, "--no-color", "--region", "g1:9+3");
  }

  public void testExpandReference() {
    final int[] inserts = {0, 1, 2, 3, 4};
    final String reference = "actga";
    assertEquals("a_c__t___g____a", Aview.expandReference(inserts, reference, new DisplayHelper()));
  }

  public void testInputList() throws IOException {
    final File f = FileUtils.createTempDir("aview", "inputlist");
    try {
      prepareData(f, SAM);
      final String path = new File(f, ALIGNMENTS_BAM_FILE_NAME).getPath();
      final File list = new File(f, "list");
      FileUtils.stringToFile(path + LS, list);
      final String[] args = {
          "-I", list.getPath(),
          "-c", new File(f, SNPS_FILE_NAME).getPath(),
          "-b", new File(f, SNPS_FILE_NAME).getPath(),
          "--region", "g1:11+3",
          "-t", new File(f, TEMPLATE_FILE_NAME).getPath(),
          "--no-dots",
          "--sort-reads",
      };
      final Aview aview = new Aview();
      assertEquals("aview", aview.moduleName());
      final ByteArrayOutputStream bos = new ByteArrayOutputStream();
      final MemoryPrintStream err = new MemoryPrintStream();
      final int code2 = aview.mainInit(args, bos, err.printStream());
      GlobalFlags.resetAccessedStatus();
      assertEquals(err.toString(), 0, code2);
      assertEquals("", err.toString());
      bos.flush();
      final String res = bos.toString();
      //System.err.println(res);
      mNano.check("aview-testinputlist.txt", res);
      bos.close();
      err.close();
    } finally {
      assertTrue(FileHelper.deleteAll(f));
    }
  }

  static final String SAM_RG = "" + "@HD" + T + "VN:1.0" + T + "SO:coordinate" + LS
  + "@SQ" + T + "SN:g1" + T + "LN:20" + LS
  + "@RG" + T + "ID:A" + T + "PLATFORM:ILLUMINA" + T + "SM:NA12" + LS
  + "@RG" + T + "ID:B" + T + "PLATFORM:ILLUMINA" + T + "SM:NA12" + LS
  + "0" + T + "0" + T + "g1" + T + "3" + T + "255" + T + "8M" + T + "*" + T + "0" + T + "0" + T + "ATCGACTG" + T + "&'(`````" + T + "AS:i:0" + T + "IH:i:1" + T + "RG:Z:B" + LS
  + "1" + T + "0" + T + "g1" + T + "3" + T + "255" + T + "8M" + T + "*" + T + "0" + T + "0" + T + "ATCGACTG" + T + "````````" + T + "AS:i:0" + T + "IH:i:1" + T + "RG:Z:B" + LS
  + "2" + T + "0" + T + "g1" + T + "5" + T + "255" + T + "8M" + T + "*" + T + "0" + T + "0" + T + "CGACTGGT" + T + "````````" + T + "AS:i:1" + T + "IH:i:1" + T + "RG:Z:B" + LS
  + "3" + T + "0" + T + "g1" + T + "6" + T + "255" + T + "8M" + T + "*" + T + "0" + T + "0" + T + "GACTGGTC" + T + "````````" + T + "AS:i:1" + T + "IH:i:1" + T + "RG:Z:A" + LS
  + "4" + T + "0" + T + "g1" + T + "6" + T + "255" + T + "8M" + T + "*" + T + "0" + T + "0" + T + "GACTGGTC" + T + "````````" + T + "AS:i:1" + T + "IH:i:1" + T + "RG:Z:A" + LS
  + "5" + T + "0" + T + "g1" + T + "11" + T + "255" + T + "8M" + T + "*" + T + "0" + T + "0" + T + "GTCAGCTA" + T + "````````" + T + "AS:i:1" + T + "IH:i:1" + T + "RG:Z:A" + LS
  ;

  public void testRG() throws IOException, InterruptedException {
    run("aview-rg.txt", SAM_RG, "--region", "g1:1+11", "--no-dots", "--sort-readgroup", "--no-color", "--print-readgroup", "--print-mapq");
  }

  public void testTooMuchTemplate() throws Exception {
    runWithErr("aview-toomuchtemplate.txt", SAM, "The end position \"23\" is outside the length of the sequence (20). Defaulting end to \"20\"" + LS, "--region", "g1:11+13", "--no-dots", "--no-base-colors");
  }
}
