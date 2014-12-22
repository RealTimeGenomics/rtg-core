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

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.tabix.TabixIndexer;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.cli.CFlags;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.test.FileHelper;
import com.rtg.vcf.SnpIntersection.LineHolder;

/**
 * Test class for {@link com.rtg.vcf.SnpIntersection}
 */
public final class SnpIntersectionTest extends AbstractCliTest {

  private static final String HEADER = ""
      + "##fileformat=VCFv4.1\n"
      + "##fileDate=20110727"
      + "##source=RTGvPOST-2.2-DEV build <not found> (<not found>)\n"
      + "##CL=null\n"
      + "##RUN-ID=74c45566-c354-42e3-b774-b84badbb604c\n"
      + "##TEMPLATE-SDF-ID=76c1612e-bb6d-4244-a5e7-6f2415c980ec\n"
      + "##FILTER=<ID=c-1,Description=\"Coverage exceeded -1\">\n"
      + "##FILTER=<ID=a1000.0,Description=\"Ambiguity exceeded 1000.0\">\n"
      + "##FILTER=<ID=RC,Description=\"RTG variant is a complex region\">\n"
      + "##FILTER=<ID=RX,Description=\"RTG variant contained within hypercomplex region\">\n"
      + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
      + "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n"
      + "##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"Genotype Quality\">\n"
      + "##FORMAT=<ID=RE,Number=1,Type=Float,Description=\"RTG Total Error\">\n"
      + "##FORMAT=<ID=AR,Number=1,Type=Float,Description=\"Ambiguity Ratio\">\n"
      + "##FORMAT=<ID=AB,Number=1,Type=Float,Description=\"Allele Balance\">\n"
      + "##FORMAT=<ID=RS,Number=.,Type=String,Description=\"RTG Support Statistics\">\n"
      + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n";
  private static final String SAME = "chr1\t147\t.\tC\tA\t.\tPASS\t%s\tGT:DP:RE:GQ:RS\t1/0:16:0.670:30.0:A,6,0.513,C,10,0.157";
  private static final String DIFFERENT1 = "chr1\t190\t.\tC\tT,G\t.\tPASS\t%s\tGT:DP:RE:GQ:RS\t1/2:16:0.670:30.0:T,6,0.513,G,10,0.157";
  private static final String DIFFERENT2 = "chr1\t190\t.\tC\tG\t.\tPASS\t%s\tGT:DP:RE:GQ:RS\t1/1:16:0.670:30.0:T,6,0.513,G,10,0.157";
  private static final String FIRST_ONLY1 = "chr1\t180\t.\tC\tT,G\t.\tPASS\t%s\tGT:DP:RE:GQ:RS\t1/2:16:0.670:30.0:T,6,0.513,G,10,0.157";
  private static final String FIRST_ONLY2 = "chr1\t255\t.\tC\tT,G\t.\tPASS\t%s\tGT:DP:RE:GQ:RS\t1/2:16:0.670:30.0:T,6,0.513,G,10,0.157";
  private static final String SECOND_ONLY = "chr1\t260\t.\tC\tT,G\t.\tPASS\t%s\tGT:DP:RE:GQ:RS\t1/2:16:0.670:30.0:T,6,0.513,G,10,0.157";

  private static final String INPUT_STR_FIRST = HEADER
      + String.format(SAME, ".") + StringUtils.LS
      + String.format(FIRST_ONLY1, ".") + StringUtils.LS
      + String.format(DIFFERENT1, ".") + StringUtils.LS
      + String.format(FIRST_ONLY2, ".") + StringUtils.LS;

  private static final String INPUT_STR_SECOND = HEADER.replaceAll("Float", "Integer") //test header merging

      + String.format(SAME, ".") + StringUtils.LS
      + String.format(DIFFERENT2, ".") + StringUtils.LS
      + String.format(SECOND_ONLY, ".") + StringUtils.LS;

  private static final String EXPECTED_MESSAGE1 =
      "Total results in first file           : 4" + StringUtils.LS
      + "Total results in second file          : 3" + StringUtils.LS
      + "Number of results same                : 1" + StringUtils.LS
      + "Number of results different           : 1" + StringUtils.LS
      + "Number of results only in first file  : 2" + StringUtils.LS
      + "Number of results only in second file : 1" + StringUtils.LS;

  private static final String EXPECTED_MESSAGE2 =
      "Total results in first file           : 3" + StringUtils.LS
      + "Total results in second file          : 4" + StringUtils.LS
      + "Number of results same                : 1" + StringUtils.LS
      + "Number of results different           : 1" + StringUtils.LS
      + "Number of results only in first file  : 1" + StringUtils.LS
      + "Number of results only in second file : 2" + StringUtils.LS;

  @Override
  protected AbstractCli getCli() {
    return new SnpIntersection();
  }

  public void testHelp() {
    checkHelp("rtg snpintersect"
        , "Produces intersection information between two SNP files."
        , "-i,", "input-first=FILE", "first file"
        , "-I,", "input-second=FILE", "second file"
        , "-o,", "output=DIR", "output directory"
        , "--region=STRING", "if set, only process the SNPs within the specified range. The format is one of <template_name>, <template_name>:start-end or <template_name>:start+length"
        , "-Z,", "no-gzip", "do not gzip the output"
        , "-c,", "compare-alts", "do a basic comparison on ALT calls in addition to position"
        );
  }

  public void testValidator() throws IOException {
    final File tempDir = FileUtils.createTempDir("snpIntersectionValidator", "test");
    try {
      final MemoryPrintStream ps = new MemoryPrintStream();
      Diagnostic.setLogStream(ps.printStream());
      final CFlags flags = new CFlags("blah", ps.printStream(), ps.printStream());
      SnpIntersection.initFlags(flags);
      assertFalse(flags.setFlags("-o", tempDir.getPath(), "-i", "blah", "-I", "blah"));
      ps.reset();
      assertFalse(flags.setFlags("-o", new File(tempDir.getPath(), "output").getPath(), "-i", "blah", "-I", "blah"));
      assertTrue(ps.toString(), ps.toString().contains("The file \"blah\", specified for \"--input-first\", does not exist."));
      ps.reset();
      final File input = new File(tempDir, "file");
      assertTrue(input.createNewFile());
      assertFalse(flags.setFlags("-o", new File(tempDir.getPath(), "output").getPath(), "-i", input.getPath(), "-I", tempDir.getPath()));
      assertTrue(ps.toString(), ps.toString().contains("The file \"" + tempDir.getPath() + "\", specified for \"--input-second\", is not a file."));
      ps.reset();
      assertFalse(flags.setFlags("-o", new File(tempDir.getPath(), "output").getPath(), "-i", input.getPath(), "-I", input.getPath(), "--region", "chr1:0-10"));
      assertTrue(ps.toString(), ps.toString().contains("The value \"chr1:0-10\" for \"--region\" is not a well formed region."));
      try {
        flags.setFlags("-o", new File(tempDir.getPath(), "output").getPath(), "-i", input.getPath(), "-I", input.getPath(), "--region", "chr1:1-10");
        fail();
      } catch (final NoTalkbackSlimException e) {
        assertEquals("\"" + input.getPath() + "\" is not in bgzip format. Cannot use \"--region\".", e.getMessage());
      }
      final File bgzipInput = new File(tempDir, "file.gz");
      FileHelper.stringToGzFile("#Test String", bgzipInput);
      try {
        flags.setFlags("-o", new File(tempDir.getPath(), "output").getPath(), "-i", bgzipInput.getPath(), "-I", input.getPath(), "--region", "chr1:1-10");
        fail();
      } catch (final NoTalkbackSlimException e) {
        assertEquals("\"" + input.getPath() + "\" is not in bgzip format. Cannot use \"--region\".", e.getMessage());
      }
      try {
        flags.setFlags("-o", new File(tempDir.getPath(), "output").getPath(), "-i", bgzipInput.getPath(), "-I", bgzipInput.getPath(), "--region", "chr1:1-10");
        fail();
      } catch (final NoTalkbackSlimException e) {
        assertEquals("Index not found for file: \"" + bgzipInput.getPath() + "\" expected index called: \"" + bgzipInput.getPath() + ".tbi\"", e.getMessage());
      }
      final File bgzipInputWithIndex = new File(tempDir, "file2.gz");
      FileHelper.stringToGzFile("#Test String", bgzipInputWithIndex);
      final File bgzipInputIndex = new File(tempDir, "file2.gz.tbi");
      assertTrue(bgzipInputIndex.createNewFile());
      try {
        flags.setFlags("-o", new File(tempDir.getPath(), "output").getPath(), "-i", bgzipInputWithIndex.getPath(), "-I", bgzipInput.getPath(), "--region", "chr1:1-10");
        fail();
      } catch (final NoTalkbackSlimException e) {
        assertEquals("Index not found for file: \"" + bgzipInput.getPath() + "\" expected index called: \"" + bgzipInput.getPath() + ".tbi\"", e.getMessage());
      }
    } finally {
      assertTrue(FileHelper.deleteAll(tempDir));
    }
  }

  public void testComparePredictions() {
    assertTrue(SnpIntersection.comparePredictions("bl:ah", "ah:bl"));
    assertTrue(SnpIntersection.comparePredictions("blah", "blah"));
    assertEquals("ah:bl", SnpIntersection.reverseDiploidPrediction("bl:ah"));
    assertEquals("blah", SnpIntersection.reverseDiploidPrediction("blah"));
  }

  public void testRegion() throws Exception {
    final File topLevel = FileUtils.createTempDir("snpintersection", "regiontest");
    try {
      final File first = new File(topLevel, "first.txt.gz");
      FileHelper.stringToGzFile(INPUT_STR_FIRST, first);
      new TabixIndexer(first, new File(topLevel, "first.txt.gz.tbi")).saveTsvIndex();
      final File second = new File(topLevel, "second.txt.gz");
      FileHelper.stringToGzFile(INPUT_STR_SECOND, second);
      new TabixIndexer(second, new File(topLevel, "second.txt.gz.tbi")).saveTsvIndex();
      final File outDir = new File(topLevel, "out");
      final String[] args = {
          "-i", first.getPath(),
          "-I", second.getPath(),
          "-o", outDir.getPath(),
          "-c",
          "--region", "chr1:148-254"
      };
      final String output = checkMainInitOk(args);
      TestUtils.containsAll(output
          , "Total results in first file           : 3"
          , "Total results in second file          : 2"
          , "Number of results same                : 1"
          , "Number of results different           : 1"
          , "Number of results only in first file  : 1"
          , "Number of results only in second file : 0"
          );
      TestUtils.containsAll(FileHelper.gzFileToString(new File(outDir, "first-only.vcf.gz")), String.format(FIRST_ONLY1, "SF=0"));
      TestUtils.containsAll(FileHelper.gzFileToString(new File(outDir, "same.vcf.gz")), String.format(SAME, "SF=0"), String.format(SAME, "SF=1"));
      TestUtils.containsAll(FileHelper.gzFileToString(new File(outDir, "different.vcf.gz")), String.format(DIFFERENT1, "SF=0"), String.format(DIFFERENT2, "SF=1"));
    } finally {
      assertTrue(FileHelper.deleteAll(topLevel));
    }
  }

  public void testMainExec() throws IOException {
    final File topLevel = FileUtils.createTempDir("snpintersection", "test");
    try {
      final File first = new File(topLevel, "first");
      FileUtils.stringToFile(INPUT_STR_FIRST, first);

      final File second = new File(topLevel, "second");
      FileUtils.stringToFile(INPUT_STR_SECOND, second);

      run(topLevel, first, second, false);
    } finally {
      assertTrue(FileHelper.deleteAll(topLevel));
    }
  }

  public void testMainExecGzipped() throws IOException {
    final File topLevel = FileUtils.createTempDir("snpintersection", "test");
    try {
      final File first = new File(topLevel, "first.gz");
      FileHelper.stringToGzFile(INPUT_STR_FIRST, first);

      final File second = new File(topLevel, "second.gz");
      FileHelper.stringToGzFile(INPUT_STR_SECOND, second);

      run(topLevel, first, second, true);
    } finally {
      assertTrue(FileHelper.deleteAll(topLevel));
    }
  }

  private void run(final File topLevel, final File first, final File second, boolean gzip) throws IOException {
    final File outDirFirstSecond = new File(topLevel, "out-forward");
    final File outDirSecondFirst = new File(topLevel, "out-reverse");

    final String gzipSuff = gzip ? FileUtils.GZ_SUFFIX : "";
    final SnpIntersection ins = new SnpIntersection();

    final String[] argsForward;
    final String[] argsReverse;
    if (gzip) {
      argsForward = new String[] {"-i", first.getPath(),
          "-I", second.getPath(),
          "-o", outDirFirstSecond.getPath(), "-c"};
      argsReverse = new String[] {"-i", second.getPath(),
          "-I", first.getPath(),
          "-o", outDirSecondFirst.getPath(), "-c"};
    } else {
      argsForward = new String[] {"-i", first.getPath(),
          "-I", second.getPath(),
          "-o", outDirFirstSecond.getPath(),
          "-Z", "-c"};
      argsReverse = new String[] {"-i", second.getPath(),
          "-I", first.getPath(),
          "-o", outDirSecondFirst.getPath(),
          "-Z", "-c"};
    }
    final ByteArrayOutputStream bos = new ByteArrayOutputStream();
    try {
      final MemoryPrintStream mps = new MemoryPrintStream();
      final int code = ins.mainInit(argsForward, bos, mps.printStream());
      assertEquals(mps.toString(), 0, code);
    } finally {
      bos.close();
    }

    //System.err.println(bos.toString());
    final String status = FileUtils.fileToString(new File(outDirFirstSecond, "summary.txt"));
    assertEquals(EXPECTED_MESSAGE1, status);
    assertEquals(EXPECTED_MESSAGE1, bos.toString());
    final File firstOnlyForward = new File(outDirFirstSecond, "first-only.vcf" + gzipSuff);
    assertTrue(firstOnlyForward.exists());
    final File secondOnlyForward = new File(outDirFirstSecond, "second-only.vcf" + gzipSuff);
    assertTrue(secondOnlyForward.exists());
    final File sameForward = new File(outDirFirstSecond, "same.vcf" + gzipSuff);
    assertTrue(sameForward.exists());
    final File differentForward = new File(outDirFirstSecond, "different.vcf" + gzipSuff);
    assertTrue(differentForward.exists());
    if (gzip) {
      assertTrue(new File(firstOnlyForward.getPath() + TabixIndexer.TABIX_EXTENSION).exists());
      assertTrue(new File(secondOnlyForward.getPath() + TabixIndexer.TABIX_EXTENSION).exists());
      assertTrue(new File(sameForward.getPath() + TabixIndexer.TABIX_EXTENSION).exists());
      assertTrue(new File(differentForward.getPath() + TabixIndexer.TABIX_EXTENSION).exists());
    }
    TestUtils.containsAll(gzip ? FileHelper.gzFileToString(firstOnlyForward) : FileUtils.fileToString(firstOnlyForward), String.format(FIRST_ONLY1, "SF=0"), String.format(FIRST_ONLY2, "SF=0"));
    TestUtils.containsAll(gzip ? FileHelper.gzFileToString(secondOnlyForward) : FileUtils.fileToString(secondOnlyForward), String.format(SECOND_ONLY, "SF=1"));
    TestUtils.containsAll(gzip ? FileHelper.gzFileToString(sameForward) : FileUtils.fileToString(sameForward), String.format(SAME, "SF=0"), String.format(SAME, "SF=1"));
    TestUtils.containsAll(gzip ? FileHelper.gzFileToString(differentForward) : FileUtils.fileToString(differentForward), String.format(DIFFERENT1, "SF=0"), String.format(DIFFERENT2, "SF=1"));
    TestUtils.containsAll(FileUtils.fileToString(new File(outDirFirstSecond, "progress")),
        "Start loading first file...",
        "Start loading second file...");

    final ByteArrayOutputStream bos2 = new ByteArrayOutputStream();
    try {
      assertEquals(0, ins.mainInit(argsReverse, bos2, TestUtils.getNullPrintStream()));
    } finally {
      bos2.close();
    }
    final String status2 = FileUtils.fileToString(new File(outDirSecondFirst, "summary.txt"));
    assertEquals(EXPECTED_MESSAGE2, status2);
    assertEquals(EXPECTED_MESSAGE2, bos2.toString());
    final File firstOnlyReverse = new File(outDirSecondFirst, "first-only.vcf" + gzipSuff);
    assertTrue(firstOnlyReverse.exists());
    final File secondOnlyReverse = new File(outDirSecondFirst, "second-only.vcf" + gzipSuff);
    assertTrue(secondOnlyReverse.exists());
    final File sameReverse = new File(outDirSecondFirst, "same.vcf" + gzipSuff);
    assertTrue(sameReverse.exists());
    final File differentReverse = new File(outDirSecondFirst, "different.vcf" + gzipSuff);
    assertTrue(differentReverse.exists());


    TestUtils.containsAll(gzip ? FileHelper.gzFileToString(firstOnlyReverse) : FileUtils.fileToString(firstOnlyReverse), String.format(SECOND_ONLY, "SF=0"));
    TestUtils.containsAll(gzip ? FileHelper.gzFileToString(secondOnlyReverse) : FileUtils.fileToString(secondOnlyReverse), String.format(FIRST_ONLY1, "SF=1"), String.format(FIRST_ONLY2, "SF=1"));
    TestUtils.containsAll(gzip ? FileHelper.gzFileToString(sameReverse) : FileUtils.fileToString(sameReverse), String.format(SAME, "SF=0"), String.format(SAME, "SF=1"));
    TestUtils.containsAll(gzip ? FileHelper.gzFileToString(differentReverse) : FileUtils.fileToString(differentReverse), String.format(DIFFERENT1, "SF=1"), String.format(DIFFERENT2, "SF=0"));
    TestUtils.containsAll(FileUtils.fileToString(new File(outDirSecondFirst, "progress")),
        "Start loading first file...",
        "Start loading second file...");
    if (gzip) {
      assertTrue(new File(firstOnlyReverse.getPath() + TabixIndexer.TABIX_EXTENSION).exists());
      assertTrue(new File(secondOnlyReverse.getPath() + TabixIndexer.TABIX_EXTENSION).exists());
      assertTrue(new File(sameReverse.getPath() + TabixIndexer.TABIX_EXTENSION).exists());
      assertTrue(new File(differentReverse.getPath() + TabixIndexer.TABIX_EXTENSION).exists());
    }
    TestUtils.containsAll(FileUtils.fileToString(new File(outDirSecondFirst, "snpintersect.log")), "progress=" + false, gzip ? "SnpIndex" : "");
  }

  public void testComparator() {
    final LineHolder holder1 = new LineHolder(null);
    final LineHolder holder2 = new LineHolder(null);
    final LineHolder holder3 = new LineHolder(null);
    final LineHolder holder4 = new LineHolder(null);
    final List<LineHolder> list = new ArrayList<>();
    list.add(holder4);
    list.add(holder2);
    list.add(holder1);
    list.add(holder3);
    Collections.sort(list);
    assertTrue(holder1 == list.get(0));
    assertTrue(holder2 == list.get(1));
    assertTrue(holder3 == list.get(2));
    assertTrue(holder4 == list.get(3));
    assertFalse(holder1.equals(false));
    assertTrue(holder1.equals(holder1));
    assertFalse(holder1.equals(holder2));
    assertFalse(holder1.equals(list));
    assertTrue(holder1.hashCode() == holder1.mInstance);
  }
}
