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
package com.rtg.protein;

import static com.rtg.util.StringUtils.LS;
import static com.rtg.util.StringUtils.TAB;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

import com.rtg.mode.DnaUtils;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 * Test for corresponding class
 */
public class MapXFunctionalTest extends TestCase {

  private static final int HEADER_LINES = 6;

  private static final String[] READS_PERFECT = {
    "AAATGGCGCAAAAACAGAAAGTCGAAAAAAAATCAA",
    "AATGGCGCAAAAACAGAAAGTCGAAAAAAAATCAAA",
    "ATGGCGCAAAAACAGAAAGTCGAAAAAAAATCAAAG",
    "TGGCGCAAAAACAGAAAGTCGAAAAAAAATCAAAGA",
    DnaUtils.reverseComplement("AAATGGCGCAAAAACAGAAAGTCGAAAAAAAATCAA"),
    DnaUtils.reverseComplement("AATGGCGCAAAAACAGAAAGTCGAAAAAAAATCAAA"),
    DnaUtils.reverseComplement("ATGGCGCAAAAACAGAAAGTCGAAAAAAAATCAAAG"),
    DnaUtils.reverseComplement("TGGCGCAAAAACAGAAAGTCGAAAAAAAATCAAAGA"),
  "AAAAAAATCAAAGAAATTATAACCACGACGCAGCAG"};
  private static final String[] READS_ONE_SUB = {
    "AAATGGCGCAAAAACAGACAGTCGAAAAAAAATCAA",
    "AATGGCGCAAAGACAGAAAGTCGAAAAAAAATCAAA",
    "ATGGCGCAAAAACAGAAAGTCGAAATAAAATCAAAG",
    "TGGCGCAAAAACAGAAAGTCGACAAAAAATCAAAGA",
    DnaUtils.reverseComplement("AAATGGCGCAAAAACAGAAAGTCGAAAAAAAATCAA"),
    DnaUtils.reverseComplement("AATGCCGCAAAAACAGAAAGTCGAAAAAAAATCAAA"),
    DnaUtils.reverseComplement("ATGGCCCAAAAACAGAAAGTCGAAAAAAAATCAAAG"),
    DnaUtils.reverseComplement("TGGCGCAATAACAGAAAGTCGAAAAAAAATCAAAGA"),
    "AAAAAAATCAAAGAAATTATAACCACGACGCGGCAG"
  };

  private static final String READS_FASTA_PERFECT;
  private static final long LENGTH_FASTA_PERFECT;
  private static final String READS_FASTA_ONE_SUB;
  private static final long LENGTH_FASTA_ONE_SUB;
  private static final String TEMPLATE_FASTA_2;
  static {
    StringBuilder sb = new StringBuilder();
    long t0 = 0;
    for (int i = 0; i < READS_PERFECT.length; i++) {
      sb.append(">testRead").append(i).append(LS).append(READS_PERFECT[i]).append(LS);
      t0 += READS_PERFECT[i].length();
    }
    READS_FASTA_PERFECT = sb.toString();
    LENGTH_FASTA_PERFECT = t0;

    sb = new StringBuilder();
    long t1 = 0;
    for (int i = 0; i < READS_ONE_SUB.length; i++) {
      sb.append(">testRead").append(i).append(LS).append(READS_ONE_SUB[i]).append(LS);
      t1 += READS_ONE_SUB[i].length();
    }
    READS_FASTA_ONE_SUB = sb.toString();
    LENGTH_FASTA_ONE_SUB = t1;

    sb = new StringBuilder();
    for (int i = 0; i < 23; i++) {
      sb.append(">templateName").append(i).append(LS).append(MapXCliTest.TEMPLATE_PROTEIN).append(LS);
    }
    TEMPLATE_FASTA_2 = sb.toString();
  }

  public void testPerfectMatch() throws IOException {
    final String expected = FileHelper.resourceToString("com/rtg/protein/resources/mapxtestresults.txt");
    check(MapXCliTest.TEMPLATE_FASTA, READS_FASTA_PERFECT, expected, new String[] {"-a", "1", "-b", "0", "-w", "9", "-T", "1"}, true, LENGTH_FASTA_PERFECT);
  }

  public void testPerfectMatchPartial() throws IOException {
    final String expected = FileHelper.resourceToString("com/rtg/protein/resources/mapxtestresults_partial.txt");
    check(MapXCliTest.TEMPLATE_FASTA, READS_FASTA_PERFECT, expected, new String[] {"-a", "1", "-b", "0", "-w", "9", "-T", "1", "--start-read", "2", "--end-read", "5"}, true, 108L);
  }

  public void testPerfectMatchMC() throws IOException {
    final File dir = FileUtils.createTempDir("mapx", "test");
    try {
      final File template = new File(dir, "template");
      final File reads = new File(dir, "reads");
      final File output = new File(dir, "output");
      final File output2 = new File(dir, "output2");
      ReaderTestUtils.getReaderProtein(TEMPLATE_FASTA_2, template).close();
      ReaderTestUtils.getReaderDNA(READS_FASTA_PERFECT, reads, null).close();
      final MapXCli foo = new MapXCli();
      final MemoryPrintStream ps = new MemoryPrintStream();
      final int rc = foo.mainInit(new String[] {"-Z", "-t", template.getPath(), "-i", reads.getPath(), "-a", "1", "-b", "0", "-o", output.getPath(), "-w", "9", "-T", "10"}, TestUtils.getNullOutputStream(), ps.printStream());
      assertEquals(ps.toString(), 0, rc);
      ps.close();
      final MapXCli foo2 = new MapXCli();
      assertEquals(0, foo2.mainInit(new String[] {"-Z", "-t", template.getPath(), "-i", reads.getPath(), "-a", "1", "-b", "0", "-o", output2.getPath(), "-w", "9", "-T", "1"}, TestUtils.getNullOutputStream(), TestUtils.getNullPrintStream()));
      final File leftFile = new File(output, "alignments.tsv");
      final File rightFile = new File(output2, "alignments.tsv");
      assertEquals(leftFile.length(), rightFile.length());
      assertTrue(rightFile.length() > 0);
      final String left = FileUtils.fileToString(leftFile);
      final String right = FileUtils.fileToString(rightFile);

      assertEquals(left, right);
      final String log = FileUtils.fileToString(new File(output, "mapx.log"));
      TestUtils.containsAll(log, "-T 10");
      //assertEquals(FileUtils.fileToString(new File(output2, "out")), FileUtils.fileToString(new File(output, "out")));
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }
  public void testOneSubMatch() throws IOException {
    final String expected = FileHelper.resourceToString("com/rtg/protein/resources/mapxtestresults_one_sub.txt");
    check(MapXCliTest.TEMPLATE_FASTA, READS_FASTA_ONE_SUB, expected, new String[] {"-a", "1", "-b", "0", "-w", "9", "-T", "1"}, true, LENGTH_FASTA_ONE_SUB);
  }

  public void testOneIndelMatch() throws IOException {
    final File dir = FileUtils.createTempDir("mapx", "test");
    try {
      final File template = new File(dir, "template");
      final File reads = new File(dir, "reads");
      final File output = new File(dir, "output");
      ReaderTestUtils.getReaderProtein(MapXCliTest.TEMPLATE_FASTA, template).close();
      ReaderTestUtils.getReaderDNA(MapXCliTest.READS_FASTA_ONE_INDEL, reads, null).close();
      final MapXCli foo = new MapXCli();
      assertEquals(0, foo.mainInit(new String[] {"-t", template.getPath(), "-i", reads.getPath(), "-a", "1", "-b", "1", "-o", output.getPath(), "-w", "4", "-T", "1"}, TestUtils.getNullOutputStream(), TestUtils.getNullPrintStream()));
      //Printout for convenience for the following error, the second print will produce the read in question in the +1 frame
      //System.out.println(TEMPLATE_PROTEIN);
      //System.out.println(SimpleTestUtils.dnaToProtein(READS_ONE_INDEL[0].substring(0)));
      TestUtils.containsAll(FileHelper.gzFileToString(new File(output, "alignments.tsv.gz")), "templateName\t+1\t0\t2");
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  public void testAlignmentScoreFiltering() throws IOException {
    final File dir = FileUtils.createTempDir("mapx", "test");
    try {
      final File template = new File(dir, "template");
      final File reads = new File(dir, "reads");
      final File output = new File(dir, "output");
      ReaderTestUtils.getReaderProtein(MapXCliTest.TEMPLATE_FASTA, template).close();
      ReaderTestUtils.getReaderDNA(MapXCliTest.READS_FASTA_MULTI_SUB, reads, null).close();
      final MapXCli foo = new MapXCli();
      CommandLine.setCommandArgs("automonkey", "angry");
      try {
        assertEquals(0, foo.mainInit(new String[] {"-Z", "-t", template.getPath(), "-i", reads.getPath(), "-a", "1", "-b", "0", "-o", output.getPath(), "-w", "9", "-T", "1", "-e", "-60"}, TestUtils.getNullOutputStream(), TestUtils.getNullPrintStream()));
        String results = FileUtils.fileToString(new File(output, "alignments.tsv"));
        assertEquals(results, HEADER_LINES + 4, results.split(LS).length);
        assertTrue(results.contains("#CL\tautomonkey angry"));
        assertTrue(FileHelper.deleteAll(output));
        assertEquals(0, foo.mainInit(new String[] {"-Z", "-t", template.getPath(), "-i", reads.getPath(), "-a", "1", "-b", "0", "-o", output.getPath(), "-w", "9", "-T", "1", "-e", "-61"}, TestUtils.getNullOutputStream(), TestUtils.getNullPrintStream()));
        results = FileUtils.fileToString(new File(output, "alignments.tsv"));
        assertEquals(HEADER_LINES + 4, results.split(LS).length);
        assertTrue(results.contains("#CL\tautomonkey angry"));

        assertTrue(FileHelper.deleteAll(output));
        assertEquals(0, foo.mainInit(new String[] {"-Z", "-t", template.getPath(), "-i", reads.getPath(), "-a", "1", "-b", "0", "-o", output.getPath(), "-w", "9", "-T", "1", "-e", "-62"}, TestUtils.getNullOutputStream(), TestUtils.getNullPrintStream()));
        results = FileUtils.fileToString(new File(output, "alignments.tsv"));
        assertTrue(results.contains("#CL\tautomonkey angry"));
        assertEquals(HEADER_LINES + 3, results.split(LS).length);

        assertTrue(FileHelper.deleteAll(output));
        assertEquals(0, foo.mainInit(new String[] {"-t", template.getPath(), "-i", reads.getPath(), "-a", "1", "-b", "0", "-o", output.getPath(), "-w", "9", "-T", "2", "-e", "-62", "-Z", "--no-unmapped"}, TestUtils.getNullOutputStream(), TestUtils.getNullPrintStream()));
        results = FileUtils.fileToString(new File(output, "alignments.tsv"));
        assertTrue(results.contains("#CL\tautomonkey angry"));
        assertFalse(new File(output, "unmapped.tsv").exists());
        assertEquals(HEADER_LINES + 3, results.split(LS).length);
      } finally {
        CommandLine.clearCommandArgs();
      }
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  public void testPercentIdentityFiltering() throws IOException {
    final File dir = FileUtils.createTempDir("mapx", "test");
    try {
      final File template = new File(dir, "template");
      final File reads = new File(dir, "reads");
      final File output = new File(dir, "output");
      ReaderTestUtils.getReaderProtein(MapXCliTest.TEMPLATE_FASTA, template).close();
      ReaderTestUtils.getReaderDNA(MapXCliTest.READS_FASTA_MULTI_SUB, reads, null).close();
      final MapXCli foo = new MapXCli();
      assertEquals(0, foo.mainInit(new String[] {"-Z", "-t", template.getPath(), "-i", reads.getPath(), "-a", "1", "-b", "0", "-o", output.getPath(), "-w", "9", "-T", "1", "-P", "0"}, TestUtils.getNullOutputStream(), TestUtils.getNullPrintStream()));
      String results = FileUtils.fileToString(new File(output, "alignments.tsv"));
      assertEquals(HEADER_LINES + 7, results.split(LS).length);

      assertTrue(FileHelper.deleteAll(output));
      assertEquals(0, foo.mainInit(new String[] {"-Z", "-t", template.getPath(), "-i", reads.getPath(), "-a", "1", "-b", "0", "-o", output.getPath(), "-w", "9", "-T", "1"}, TestUtils.getNullOutputStream(), TestUtils.getNullPrintStream()));
      results = FileUtils.fileToString(new File(output, "alignments.tsv"));
      assertEquals(HEADER_LINES + 7, results.split(LS).length);

      assertTrue(FileHelper.deleteAll(output));
      assertEquals(0, foo.mainInit(new String[] {"-Z", "-t", template.getPath(), "-i", reads.getPath(), "-a", "1", "-b", "0", "-o", output.getPath(), "-w", "9", "-T", "1", "-P", "90"}, TestUtils.getNullOutputStream(), TestUtils.getNullPrintStream()));
      results = FileUtils.fileToString(new File(output, "alignments.tsv"));
      assertEquals(HEADER_LINES + 6, results.split(LS).length);

      assertTrue(FileHelper.deleteAll(output));
      assertEquals(0, foo.mainInit(new String[] {"-t", template.getPath(), "-i", reads.getPath(), "-a", "1", "-b", "0", "-o", output.getPath(), "-w", "9", "-T", "2", "-P", "100", "-Z"}, TestUtils.getNullOutputStream(), TestUtils.getNullPrintStream()));
      results = FileUtils.fileToString(new File(output, "alignments.tsv"));
      assertEquals(HEADER_LINES + 1, results.split(LS).length);
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  public void testEScoreFiltering() throws IOException {
    final File dir = FileUtils.createTempDir("mapx", "test");
    try {
      final File template = new File(dir, "template");
      final File reads = new File(dir, "reads");
      final File output = new File(dir, "output");
      ReaderTestUtils.getReaderProtein(MapXCliTest.TEMPLATE_FASTA, template).close();
      ReaderTestUtils.getReaderDNA(MapXCliTest.READS_FASTA_MULTI_SUB, reads, null).close();
      final MapXCli foo = new MapXCli();
      assertEquals(0, foo.mainInit(new String[] {"-Z", "-t", template.getPath(), "-i", reads.getPath(), "-a", "1", "-b", "0", "-o", output.getPath(), "-w", "9", "-T", "1", "-E", "100"}, TestUtils.getNullOutputStream(), TestUtils.getNullPrintStream()));
      String results = FileUtils.fileToString(new File(output, "alignments.tsv"));
      assertEquals(HEADER_LINES + 7, results.split(LS).length);

      assertTrue(FileHelper.deleteAll(output));
      assertEquals(0, foo.mainInit(new String[] {"-Z", "-t", template.getPath(), "-i", reads.getPath(), "-a", "1", "-b", "0", "-o", output.getPath(), "-w", "9", "-T", "1"}, TestUtils.getNullOutputStream(), TestUtils.getNullPrintStream()));
      results = FileUtils.fileToString(new File(output, "alignments.tsv"));
      assertEquals(HEADER_LINES + 7, results.split(LS).length);

      assertTrue(FileHelper.deleteAll(output));
      assertEquals(0, foo.mainInit(new String[] {"-Z", "-t", template.getPath(), "-i", reads.getPath(), "-a", "1", "-b", "0", "-o", output.getPath(), "-w", "9", "-T", "1", "-E", "1.0E-06"}, TestUtils.getNullOutputStream(), TestUtils.getNullPrintStream()));
      results = FileUtils.fileToString(new File(output, "alignments.tsv"));
      assertEquals(HEADER_LINES + 6, results.split(LS).length);

      assertTrue(FileHelper.deleteAll(output));
      assertEquals(0, foo.mainInit(new String[] {"-t", template.getPath(), "-i", reads.getPath(), "-a", "1", "-b", "0", "-o", output.getPath(), "-w", "9", "-T", "2", "-E", "0.000000016", "-Z"}, TestUtils.getNullOutputStream(), TestUtils.getNullPrintStream()));
      results = FileUtils.fileToString(new File(output, "alignments.tsv"));

      assertEquals(HEADER_LINES + 1, results.split(LS).length);
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  public void testUnmappedFlags() throws IOException {
    final File dir = FileUtils.createTempDir("mapx", "test");
    try {
      final File template = new File(dir, "template");
      final File reads = new File(dir, "reads");
      final File output = new File(dir, "output");
      final File output2 = new File(dir, "output2");
      ReaderTestUtils.getReaderProtein(MapXCliTest.TEMPLATE_FASTA, template).close();
      ReaderTestUtils.getReaderDNA(MapXCliTest.READS_FASTA_MULTI_SUB, reads, null).close();
      final MapXCli foo = new MapXCli();
      assertEquals(0, foo.mainInit(new String[] {"-t", template.getPath(), "-i", reads.getPath(), "-a", "1", "-b", "0", "-o", output.getPath(), "-w", "9", "-T", "2", "-E", "0.000000016", "-Z"}, TestUtils.getNullOutputStream(), TestUtils.getNullPrintStream()));
      final String results = FileUtils.fileToString(new File(output, "alignments.tsv"));
      assertEquals(HEADER_LINES + 1, results.split(LS).length);
      assertTrue(new File(output, "unmapped.tsv").exists());
      final String expected = "#read-id\treason-unmapped" + LS
          + "0" + LS + "1" + LS
          + "2\tg" + LS
          + "3\tg" + LS
          + "5\tg" + LS
          + "6\tg" + LS
          + "7\tg" + LS
          + "8\tg" + LS;
      MapXCliTest.checkUnmappedNoHeader(expected, new File(output, "unmapped.tsv"));

      assertEquals(0, foo.mainInit(new String[] {"-t", template.getPath(), "-i", reads.getPath(), "-a", "1", "-b", "0", "-o", output2.getPath(), "-w", "9", "-T", "2", "-E", "0.000000016"}, TestUtils.getNullOutputStream(), TestUtils.getNullPrintStream()));
      final String actual = FileHelper.gzFileToString(new File(output2, "unmapped.tsv.gz"));
      assertEquals(expected, actual.substring(actual.indexOf("#read-id")));
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  public void testBitScoreFiltering() throws IOException {
    final File dir = FileUtils.createTempDir("mapx", "test");
    try {
      final File template = new File(dir, "template");
      final File reads = new File(dir, "reads");
      final File output = new File(dir, "output");
      ReaderTestUtils.getReaderProtein(MapXCliTest.TEMPLATE_FASTA, template).close();
      ReaderTestUtils.getReaderDNA(MapXCliTest.READS_FASTA_MULTI_SUB, reads, null).close();
      final MapXCli foo = new MapXCli();
      assertEquals(0, foo.mainInit(new String[] {"-Z", "-t", template.getPath(), "-i", reads.getPath(), "-a", "1", "-b", "0", "-o", output.getPath(), "-w", "9", "-T", "1", "-B", "0"}, TestUtils.getNullOutputStream(), TestUtils.getNullPrintStream()));
      String results = FileUtils.fileToString(new File(output, "alignments.tsv"));
      assertEquals(HEADER_LINES + 7, results.split(LS).length);

      assertTrue(FileHelper.deleteAll(output));
      assertEquals(0, foo.mainInit(new String[] {"-Z", "-t", template.getPath(), "-i", reads.getPath(), "-a", "1", "-b", "0", "-o", output.getPath(), "-w", "9", "-T", "1"}, TestUtils.getNullOutputStream(), TestUtils.getNullPrintStream()));
      results = FileUtils.fileToString(new File(output, "alignments.tsv"));
      assertEquals(HEADER_LINES + 7, results.split(LS).length);

      assertTrue(FileHelper.deleteAll(output));
      assertEquals(0, foo.mainInit(new String[] {"-Z", "-t", template.getPath(), "-i", reads.getPath(), "-a", "1", "-b", "0", "-o", output.getPath(), "-w", "9", "-T", "1", "-B", "22"}, TestUtils.getNullOutputStream(), TestUtils.getNullPrintStream()));
      results = FileUtils.fileToString(new File(output, "alignments.tsv"));
      assertEquals(HEADER_LINES + 6, results.split(LS).length);

      assertTrue(FileHelper.deleteAll(output));
      assertEquals(0, foo.mainInit(new String[] {"-t", template.getPath(), "-i", reads.getPath(), "-a", "1", "-b", "0", "-o", output.getPath(), "-w", "9", "-T", "2", "-B", "30.3", "-Z"}, TestUtils.getNullOutputStream(), TestUtils.getNullPrintStream()));
      results = FileUtils.fileToString(new File(output, "alignments.tsv"));
      assertEquals(HEADER_LINES + 1, results.split(LS).length);
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  private static final String HEADER = "#template-name\tframe\tread-id\ttemplate-start\ttemplate-end\ttemplate-length\tread-start\tread-end\tread-length\ttemplate-protein\tread-protein\talignment\tidentical\t%identical\tpositive\t%positive\tmismatches\traw-score\tbit-score\te-score" + LS;
  private static final String BUG_R = "AATTAAATATTACCGGTGGATGCAACGTACAGTACGCACTGCATCCGGAAACCTTCGAATATTGTGTNATCGAGG";

  private static final String BUG_DNA = ">HWUSI-EAS575_1:1:1:24:19#0/1" + LS + BUG_R + LS;
  private static final String BUG_PROT = ">gi|197303253|ref|ZP_03168294.1|" + LS
      + "GVHTGDSIVVAPSQTLGDKEYQMLRTSALNIISELNITGGCNVQYALHPETFEYCVIEVNPRVSRSSALASKATGYPIAK" + LS;

  private static final String BUG_EXPECTED = "gi|197303253|ref|ZP_03168294.1|\t+3\t0\t35\t58\t80\t3\t74\t75\tlnitggcnvqyalhpetfeycvie\tlnitggcnvqyalhpetfeycvie\tlnitggcnvqyalhpetfeycvie\t24\t100\t24\t100\t0\t-135\t56.6\t7.3e-16" + LS;
  /**
   * Now we do support 'X'
   */
  public void testBug() throws IOException {
    final String[] params = {"-a", "4", "-b", "0", "-w", "8", "-T", "1"};
    check(BUG_PROT, BUG_DNA, HEADER + BUG_EXPECTED, params, false, BUG_R.length());
  }

  private static final String BUG_R2 = "ATTCGATGAAGAGTATCCAGATGAAGATGAAGAAGACGGAGAGGCGGAGGTAGACGATGAAGACGAGTAAGATAT";
  private static final String BUG_DNA2 = ">HWUSI-EAS575_1:1:1:26:284#0/1" + LS
      + BUG_R2 + LS;

  private static final String BUG_PROT2 = ">gi|66801781|ref|XP_629810.1|" + LS
      + "MPSYNINNNQKNSKKMSSSSSSSSSSSPSSSSSSSSSSSPSSSSTTTSTTTTTSTTPTSLSPSQTSQQQQQQQQQQSTSPNLSRNNSFQFLHSFYDTFFSKLSKSGVPPSNVLQTSSTSVFSAPSSSLLTSNSLNSIGSLISPLSPTSSSVLSNSTSSIPYDGRAMNDMNVIELSPRLKALTNFKYH" + LS;

  private static final String BUG_EXPECTED2 = "gi|66801781|ref|XP_629810.1|\t-1\t0\t16\t40\t187\t1\t75\t75\tmssssssssssspsssssssssssp\tisyssssstsaspssssssgysssn\t+s sssss+s+spssssss  sss \t18\t72\t21\t84\t4\t-72\t32.3\t3.4e-8" + LS;

  public void testBug2() throws IOException  {
    final String[] params = {"-a", "4", "-b", "0", "-w", "5", "-T", "1", "-f", "topn", "-n", "2", "--Xprefilter-min-score", "25"};
    check(BUG_PROT2, BUG_DNA2, HEADER + BUG_EXPECTED2, params, false, BUG_R2.length());
  }

  private void check(final String prot, final String dna, final String expected, final String[] params, final boolean zipped, final long usageMetric) throws IOException {
    final File dir = FileUtils.createTempDir("mapx", "test");
    try {
      final File template = new File(dir, "template");
      final File reads = new File(dir, "reads");
      final File output = new File(dir, "output");
      ReaderTestUtils.getReaderProtein(prot, template).close();
      ReaderTestUtils.getReaderDNA(dna, reads, null).close();
      final MapXCli mapXCli = new MapXCli();
      final String results;
      if (zipped) {
        final String[] args = {"-t", template.getPath(), "-i", reads.getPath(), "-o", output.getPath()};
        final MemoryPrintStream pipes = new MemoryPrintStream();
        final int code = mapXCli.mainInit(TestUtils.append(args, params), pipes.outputStream(), pipes.printStream());
        assertEquals(pipes.toString(), 0, code);
        results = FileHelper.gzFileToString(new File(output, "alignments.tsv.gz"));
      } else {
        final String[] args = {"-Z", "-t", template.getPath(), "-i", reads.getPath(), "-o", output.getPath()};
        assertEquals(0, mapXCli.mainInit(TestUtils.append(args, params), TestUtils.getNullOutputStream(), TestUtils.getNullPrintStream()));
        results = FileUtils.fileToString(new File(output, "alignments.tsv"));
      }

      //System.err.println(FileUtils.fileToString(new File(output, "mapx.log")));
      if (expected != null) {
        assertTrue(TestUtils.sameLines(expected, results.substring(results.indexOf("#template-name")), true));
      }
      final String usageLog = mapXCli.usageLog();
      //System.err.println(usageLog);
      TestUtils.containsAll(usageLog, "[Usage beginning module=mapx runId=", ", Usage end module=mapx runId=", " metric=" + usageMetric + " success=true]");

    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  private static final String REV_FRAME_R = "CCTACATTTTTCTACAACCCGTGAATGGTTTTCATCAGTTGTTCTCCCAGTCCGATAGTATATCCTTGAGATACA";
  private static final String REV_FRAME_DNA = ">HWUSI-EAS575_1:1:1:10:424#0/1" + LS + REV_FRAME_R;

  private static final String REV_EXPECTED = "prot\t-2\t0\t80\t103\t94\t3\t74\t75\tvsqgytiglgeqlmkxxxxxxxxx\tvsqgytiglgeqlmktihgl*knv\tvsqgytiglgeqlmk         \t15\t63\t15\t63\t9\t-64\t29.3\t1.5e-7" + LS;
  private static final String REV_FRAME_PROT = ">prot" + LS
      + "AALGSEFEKWGRKIVLLFPSEEQYKKFHPAEFQGLPSTIAYGVDIDNRIQNQIAKSMNLSNSDILPMFIIADTFNRVVFVSQGYTIGLGEQLMK" + LS;

  public void testOffRightEndReverseFrame() throws IOException {
    final String[] params = {"-a", "5", "-b", "0", "-w", "4", "-T", "1"};
    check(REV_FRAME_PROT, REV_FRAME_DNA, HEADER + REV_EXPECTED, params, true, REV_FRAME_R.length());
  }

  private static final String REV_FRAME_R2 = "ACTACTACTACTACTACT" + "AACCCGTGAATGGTTTTCATCAGTTGTTCTCCCAGTCCGATAGTATATCCTTGAGATACA";
  private static final String REV_FRAME_DNA2 = ">HWUSI-EAS575_1:1:1:10:424#0/1" + LS + REV_FRAME_R2;

  private static final String REV_EXPECTED2 = "prot\t-2\t0\t-5\t19\t46\t3\t77\t78\txxxxxxiglgeqlmktihglaalgs\tvsqgytiglgeqlmktihglvvvvv\t      iglgeqlmktihgl  +  \t14\t56\t15\t60\t10\t-63\t28.9\t9.3e-8" + LS;
  private static final String REV_FRAME_PROT2 = ">prot" + LS
      + "IGLGEQLMKTIHGLAALGSEFEKWGRKIVLLFPSEEQYKKFHPAEF"; //QGLPSTIAYGVDIDNRIQNQIAKSMNLSNSDILPMFIIADTFNRVVF" + LS;

  public void testOffLeftEndReverseFrame() throws IOException {
    final String[] params = {"-a", "5", "-b", "0", "-w", "4", "-T", "1", "--Xprefilter-min-score", "50", "--Xprefilter-min-overlap", "50", "-P", "50"};
    check(REV_FRAME_PROT2, REV_FRAME_DNA2, HEADER + REV_EXPECTED2, params, false, REV_FRAME_R2.length());
  }

  private static final String TEMPLATE = ""
      + ">gi|83855433|ref|ZP_00948962.1|" + LS
      + "MAKETYDRSKPHLNVGTIGHVDHGKTTLTAAITKVLADAGFSEASAFDQIDNAPEEKERGITINSSHVEYATANRHYAHVDCPGHADYVKNMVTGAAQMDGAILVVAATDGPMPQTREHILLGRQVGIPRIVVFMNKVDMVDDEELLELVEMEIRDLLSFYEYDGDNGPVIQGSALGALNGEQKWVDTVLSLMEAVDSWIEEPQREVDKPFLMPIEDVFSITGRGTVATGRIETGIANTGDPVEIIGMGAEKLTSTITGIEMFRQILDRGEAGDNAGILLRGIEKSQISRGMVIVKPGSVTPHKKFKAEVYILKKEEGGRHTPFHNNYRPQFYVRTTDVTGTISLPDGVEMVMPGDNLTITVELLQTIAMNVGLRFAVREGGRTVGAGQVTEILD" + LS
      + ">gi|86131810|ref|ZP_01050407.1|" + LS
      + "MAKETYDRSKPHLNVGTIGHVDHGKTTLTAAITKVLADAGFSEASAFDQIDNAPEEKERGITINSSHVEYATENRHYAHVDCPGHADYVKNMVTGAAQMDGAILVVAATDGPMPQTREHILLGRQVGIPRIVVFMNKVDMVDDEELLELVEMEIRDLLSFYEYDGDNGPVVAGSALGALNGEQKWVDTVLELMAAVDSWIEEPLRETEKDFLMPIEDVFSITGRGTVATGRIETGIANTGDPVEIIGMGAEKLTSTITGIEMFRQILDRGEAGDNAGILLRGIEKSQISRGMVIVKPGSVTPHAKFKAEVYILKKEEGGRHTPFHNNYRPQFYVRTTDVTGNIGLPDGIEMVMPGDNLTITVELIQPIALNIGLRFAVREGGRTVGAGQVTEILD" + LS
      + ">gi|86134008|ref|ZP_01052590.1|" + LS
      + "MAKGTFDRSKPHLNIGTIGHVDHGKTTLTAAITKVLADAGFSEARSFDQIDNAPEEKERGITINTSHVEYQTANRHYAHVDCPGHADYVKNMVTGAAQMDGAILVVAATDGPMPQTREHILLGRQVGIPRIVVFLNKVDMVDDEELLELVDMEVRELLSFYEYDGDNGPVVSGSALGALNGEEKWVNTVLELMEQVDAWIEEPLREVDKDFLMPVEDVFSITGRGTVATGRIETGIANTGDVVDIIGMGAEKMSSTITGIEMFRQILDRGEAGDNAGILLRGIAKEDIKRGMVICKPGSVTPHAKFKAEVYVLKKEEGGRHTPFHNNYRPQFYVRTTDVTGTINLPSGIEMVMPGDNLTITVDLIQPIALNVGLRFAIREGGRTVGAGQVTELLD"  + LS
      + ">gi|88802460|ref|ZP_01117987.1|" + LS
      + "MAKANFDRSKPHLNIGTIGHVDHGKTTLTAAITKVLADAGFSAALSFDQIDNAPEEKERGITINTSHVEYQTANRHYAHVDCPGHADYVKNMVTGAAQMDGAILVVAATDGPMPQTREHILLGRQVGIPRMVVFMNKVDMVDDEELIELVDMEIRELLSFYEYDGDNGPVIAGSALGALNGEQKWVDTVLELMAAVDVWIEEPLREVDKPFLMPVEDVFSITGRGTVATGRIETGIANTGDTVDIIGMGAEKMTSTVTGIEMFRQILNRGEAGDNAGILLRGIAKEDIKRGMVICKPGSVTPHAKFKAEVYVLKKEEGGRHTPFHNNYRPQFYVRTTDVTGTINLPDGVEMVMPGDNLTITVDLLQPIALNIGLRFAIREGGRTVGAGQVTEVLD"  + LS
      + ">gi|260063564|ref|YP_003196644.1|" + LS
      + "MAKETFDRSKPHLNIGTIGHVDHGKTTLTAAITTVLANAGLSDIRSFDSIDNAPEEKERGITINTSHVEYQTENRHYAHVDCPGHADYVKNMVTGAAQMDGAILVVAATDGPMPQTREHILLGRQVGIPRIVVFLNKVDMVDDEELLELVEMEVRELLSFYEYDGDNSPVISGSALGALNGEQKWVDTVMELMKAVDEWIELPQRDVDKDFLMPVEDVFTITGRGTVATGRIETGVASTGDAVDIIGMGAEKLSSTITGVEMFRKILDRGEAGDNVGILLRGIEKKDIKRGMVICKPGSVTPHAKFEAEVYILKKEEGGRHTPFHNNYRPQFYVRTTDVTGTINLPDGVEMVMPGDNLTITVDLIQPIALNVGLRFAIREGGRTVGAGQVTKILD" + LS
      ;

  private static final String READ_NAME = ">HWI-EAS284_61BKE:1:1:2:1406#0/1";

  private static final String READ_BODY = "NTTATATTGGTTATTAGTCAAGAATTTCAGTAATCTGACCAGAACCTACTGTACGACCACCTTCACGGATAGCGAAACGCAAACCTACGTTAAATGCTAA";

  private static final String READ = ""
      + READ_NAME + LS
      + READ_BODY + LS
      ;

  private static final String EXPECTED_BUG1 = ""
      + "gi|86134008|ref|ZP_01052590.1|" + TAB + "-1" + TAB + "0" + TAB + "368" + TAB + "400" + TAB + "395" + TAB + "2" + TAB + "100" + TAB + "100" + TAB + "ialnvglrfaireggrtvgagqvtelldxxxxx" + TAB + "lafnvglrfaireggrtvgsgqiteild**pi*" + TAB + "+a nvglrfaireggrtvg+gq+te+ld     " + TAB + "23" + TAB + "70" + TAB + "27" + TAB + "82" + TAB + "6" + TAB + "-109" + TAB + "46.6" + TAB + "1.9e-11" + LS
      + "gi|88802460|ref|ZP_01117987.1|" + TAB + "-1" + TAB + "0" + TAB + "368" + TAB + "400" + TAB + "395" + TAB + "2" + TAB + "100" + TAB + "100" + TAB + "ialniglrfaireggrtvgagqvtevldxxxxx" + TAB + "lafnvglrfaireggrtvgsgqiteild**pi*" + TAB + "+a n+glrfaireggrtvg+gq+te+ld     " + TAB + "22" + TAB + "67" + TAB + "27" + TAB + "82" + TAB + "6" + TAB + "-109" + TAB + "46.6" + TAB + "1.9e-11" + LS
      + "gi|260063564|ref|YP_003196644.1|" + TAB + "-1" + TAB + "0" + TAB + "368" + TAB + "400" + TAB + "395" + TAB + "2" + TAB + "100" + TAB + "100" + TAB + "ialnvglrfaireggrtvgagqvtkildxxxxx" + TAB + "lafnvglrfaireggrtvgsgqiteild**pi*" + TAB + "+a nvglrfaireggrtvg+gq+t+ild     " + TAB + "23" + TAB + "70" + TAB + "27" + TAB + "82" + TAB + "6" + TAB + "-107" + TAB + "45.8" + TAB + "3.2e-11" + LS
      ;

  public void testBug1() throws IOException {
    final String[] params = {"-a", "2", "-b", "1", "-w", "8", "-T", "1"};
    check(TEMPLATE, READ, HEADER + EXPECTED_BUG1, params, true, READ_BODY.length());
  }

  public void testReductionBiggerThanHashBits() throws IOException {
    final String[] params = {"-a", "5", "-b", "1", "-w", "1", "-T", "1"};
    long length = 0;
    final StringBuilder sb = new StringBuilder();
    for (int i = 0; i < 16; i++) {
      sb.append(READ_NAME + "-").append(i).append(LS).append(READ_BODY).append(LS);
      length += READ_BODY.length();
    }
    check(TEMPLATE, sb.toString(), null, params, false, length);
  }

  private static final String READ_DNA = ">read1" + LS + "GCATCATGCACAACGTACCTTCATGCAACTTCA" + LS;

  private static final String READ_PROTEIN = ">read1+1" + LS + "ASCTTLHATS" + LS
      + ">read1+2" + LS + "HHAQRFMQL" + LS
      + ">read1+3" + LS + "IMHVPSCNF" + LS
      + ">read1-1" + LS + "*SCMKVCA*C" + LS
      + ">read1-2" + LS + "EVA*RYVHD" + LS
      + ">read1-3" + LS + "KLHEGTCMM" + LS;

  private static final String READ_PROTEIN_TOP_N = ">read1+1" + LS + "ASCTTLHATS" + LS
      + ">read1+2" + LS + "ASCTTLHATS" + LS
      + ">read1+3" + LS + "ASCTTLHATS" + LS
      + ">read1-1" + LS + "ASCTTLHATS" + LS
      + ">read1-2" + LS + "ASCTTLHATS" + LS
      + ">read1-3" + LS + "ASCTTLHATS" + LS;

  private static final String UNMAPPED_HEADER = "#read-id\treason-unmapped" + LS;

  public void testUnmapped() throws IOException {
    checkUnmapped(READ_PROTEIN, new String[] {"-e", "-60"}, UNMAPPED_HEADER + "0\td" + LS); //alignment score
    checkUnmapped(READ_PROTEIN, new String[] {"-B", "60"}, UNMAPPED_HEADER + "0\th" + LS); //bit score
    checkUnmapped(READ_PROTEIN, new String[] {"-E", "0.000000006"}, UNMAPPED_HEADER + "0\tg" + LS); //e score
    checkUnmapped(READ_PROTEIN, new String[] {"-P", "100", "--Xprefilter-min-score", "50", "--Xprefilter-min-overlap", "50"}, UNMAPPED_HEADER + "0\tf" + LS); //percent id

    //topn and tope check
    checkUnmapped(READ_PROTEIN_TOP_N, new String[] {"-f", "topn", "-n", "5"}, UNMAPPED_HEADER + "0\te" + LS);
    checkUnmapped(READ_PROTEIN_TOP_N, new String[] {"-f", "topequal", "-n", "5", "-T", "2"}, UNMAPPED_HEADER + "0\te" + LS);
  }

  private static void checkUnmapped(final String protein, final String[] params, final String expected) throws IOException {
    final File dir = FileUtils.createTempDir("mapx", "unmapped");
    try {
      final File template = new File(dir, "template");
      final File reads = new File(dir, "reads");
      final File output = new File(dir, "output");
      ReaderTestUtils.getReaderProtein(protein, template).close();
      ReaderTestUtils.getReaderDNA(READ_DNA, reads, null).close();
      final MapXCli foo = new MapXCli();
      final String[] args = {"-Z", "-t", template.getPath(), "-i", reads.getPath(), "-o", output.getPath(), "-w", "3", "-a", "1", "-b", "0"};
      final ByteArrayOutputStream errBaos = new ByteArrayOutputStream();
      final int code;
      try (PrintStream errStr = new PrintStream(errBaos)) {
        code = foo.mainInit(TestUtils.append(args, params), TestUtils.getNullOutputStream(), errStr);
      }
      assertEquals(errBaos.toString(), 0, code);
      //final String results = FileUtils.fileToString(new File(output, "alignments.tsv"));
      MapXCliTest.checkUnmappedNoHeader(expected, new File(output, "unmapped.tsv"));
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }
}
