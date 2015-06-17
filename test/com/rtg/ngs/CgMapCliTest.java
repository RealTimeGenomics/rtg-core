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
import java.io.PrintStream;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.reader.PrereadArm;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.sam.RecordIterator;
import com.rtg.sam.SkipInvalidRecordsIterator;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.test.FileHelper;

import htsjdk.samtools.SAMRecord;

/**
 * Test {@link CgMapCli}
 *
 */
public class CgMapCliTest extends AbstractCliTest {

  private static final String APP_NAME = "rtg cgmap";

  @Override
  protected AbstractCli getCli() {
    return new CgMapCli();
  }

  @Override
  public void testApplicationName() {
    assertEquals(APP_NAME, new CgMapCli().applicationName() + " " +  new CgMapCli().moduleName());
  }

  public void testCgMapFlags() {
    checkHelp("-i SDF|FILE -o DIR -t SDF",
        "reads",
        "template",
        "mask to apply",
        "maximum permitted fragment size when mating",
        "maximum mismatches allowed for mated",
        "maximum mismatches allowed for unmated",
        "minimum permitted fragment size when mating",
        "do not gzip the output",
        "maximum repeat frequency",
        "report unmated",
        "report unmapped",
        "number of threads",
        "file containing a single valid read group SAM header line"
        );
  }

  public void testCreateParams() throws Exception {
    final MemoryPrintStream err = new MemoryPrintStream();
    Diagnostic.setLogStream(err.printStream());
    final File mainOut = FileUtils.createTempDir("cgmap", "test");
    try {
      final File left = new File(mainOut, "left");
      final String inputDnaSequence = "@test\nacgtacgtacgtacgtacgtacgtacgtacgtacg\n+\n###################################";
      final String inputDnaSequence2 = "@test\nacgtacgtacgtacgtacgtacgtacgtacgtacg\n+\n###################################\n@test2\nacgtacgtacgtacgtacgtacgtacgtacgtacg\n+\n###################################";
      ReaderTestUtils.getReaderDNAFastqCG(inputDnaSequence, left, PrereadArm.LEFT).close();
      final File right = new File(mainOut, "right");
      ReaderTestUtils.getReaderDNAFastqCG(inputDnaSequence, right, PrereadArm.RIGHT).close();
      final File template = new File(mainOut, "template");
      ReaderTestUtils.getReaderDNA(">template\nacgt", template, null).close();
      err.reset();
      final File out = new File(mainOut, "out");
      final CgMapCli map = (CgMapCli) mCli;

      final String[] args = {"-i", mainOut.getPath(), "-t", template.getPath(), "-o", out.getPath(), "--end-read", "2", "--" + MapFlags.N_AS_MISMATCH};
      checkHandleFlagsOut(args);

      try (NgsParams params = map.makeParams()) {
        final String s = err.toString();
        assertTrue(s.contains("The end sequence id \"2\" is out of range, it must be from \"1\" to \"1\". Defaulting end to \"1\""));
        assertEquals(left.getPath(), params.buildFirstParams().directory().getPath());
        assertEquals(right.getPath(), params.buildSecondParams().directory().getPath());
        assertEquals(template.getPath(), params.searchParams().directory().getPath());
        assertEquals(out.getPath(), params.outputParams().directory().getPath());

        assertEquals(1000, params.maxFragmentLength().intValue());
        assertEquals(0, params.minFragmentLength().intValue());
        assertEquals(5, params.hashCountThreshold().intValue());
        assertEquals(true, params.useProportionalHashThreshold());
        assertEquals(true, params.outputParams().progress());
        assertEquals(false, params.outputParams().sorted());
        assertEquals(true, params.outputParams().isCompressOutput());
        assertEquals(false, params.outputParams().useids());
        assertEquals(false, params.outputParams().exclude());

        assertEquals(false, params.useLongReadMapping());
        assertEquals(5, params.outputParams().topN());
        assertEquals(4095, params.outputParams().errorLimit());
        assertEquals(1, params.unknownsPenalty());
        assertTrue(params.maskParams().toString().contains("Mask:cgmask"));
      }


      assertTrue(!left.exists() || FileHelper.deleteAll(left));
      ReaderTestUtils.getReaderDNA(">r0\nacgt", left, null).close();
      checkHandleFlagsOut("-i", mainOut.getPath(), "-t", template.getPath(), "-o", out.getPath());
      try {
        map.makeParams();
        fail("Expected failure on non-cg SDF");
      } catch (final InvalidParamsException e) {
        // Expected
      }
      assertTrue(!left.exists() || FileHelper.deleteAll(left));
      ReaderTestUtils.getReaderDNAFastqCG(inputDnaSequence, left, PrereadArm.LEFT);
      assertTrue(!right.exists() || FileHelper.deleteAll(right));
      ReaderTestUtils.getReaderDNA(">r0\nacgt", right, null).close();
      checkHandleFlagsOut("-i", mainOut.getPath(), "-t", template.getPath(), "-o", out.getPath());
      try {
        map.makeParams();
        fail("Expected failure on non-cg SDF");
      } catch (final InvalidParamsException e) {
        // Expected
      }
      assertTrue(!left.exists() || FileHelper.deleteAll(left));
      ReaderTestUtils.getReaderDNAFastqCG(inputDnaSequence, left, PrereadArm.LEFT).close();
      assertTrue(!right.exists() || FileHelper.deleteAll(right));
      ReaderTestUtils.getReaderDNAFastqCG(inputDnaSequence2, right, PrereadArm.RIGHT).close();
      checkHandleFlagsOut("-i", mainOut.getPath(),
          "-t", template.getPath(), "-o", out.getPath());
      try {
        err.reset();
        map.makeParams();
        fail("Expected failure on non-equal num sequences for paired-end");
      } catch (final InvalidParamsException e) {
        assertTrue(e.getMessage().contains("Left and right SDFs for read pair must have same number of sequences, actually had: 1 and 2"));
//        assertTrue(err.toString().contains("Left and right SDFs for read pair must have same number of sequences, actually had: 1 and 2"));
      }
    } finally {
      Diagnostic.setLogStream();
      assertTrue(FileHelper.deleteAll(mainOut));
      err.close();
    }
  }


  private static final String TEMPLATE = ">template" + LS
      //1234567890123456789012345678901234567890
      + "tacgtnnnnncatgactgctgcatactgcatgcatgactg"
      + "actgactgcatatgcattcatactcatgatgcatgctgca"
      + "tgatctgcatactatcattttcatgcagtactgcatgcat"
      + "gcatgnnnnncatgctactgtacgatcgatcgatcgatcg";

  private static final String R_LEFT = "tacgtcatgactgctgcatactgcatgcatgactg";
  private static final String READ_LEFT = "@left" + LS
      + R_LEFT + LS
      + "+left" + LS
      + "###################################";
  private static final long L_LENGTH = R_LEFT.length();

  private static final String R_RIGHT = "gcatgcatgctactgtacgatcgatcgatcgatcg";
  private static final String READ_RIGHT = "@right" + LS
      + R_RIGHT + LS
      + "+right" + LS
      + "###################################";
  private static final long R_LENGTH = R_RIGHT.length();
  private static final long LENGTH = L_LENGTH + R_LENGTH;


  public void testEnd2End() throws IOException {
    Diagnostic.setLogStream();
    final File outer = FileUtils.createTempDir("cgmap", "end2end");
    try {
      final File template = new File(outer, "template");
      final File reads = new File(outer, "Reads");
      final File left = new File(reads, "left");
      final File right = new File(reads, "right");
      final File out = new File(outer, "out");

      final File header = new File(outer, "header");
      FileUtils.stringToFile("@RG\tID:L23\tPL:COMPLETE\tSM:NA123", header);
      ReaderTestUtils.getReaderDNA(TEMPLATE, template, null).close();
      ReaderTestUtils.getReaderDNAFastqCG(READ_LEFT, left, PrereadArm.LEFT).close();
      ReaderTestUtils.getReaderDNAFastqCG(READ_RIGHT, right, PrereadArm.RIGHT).close();

      final CgMapCli map = new CgMapCli();

      final String[] args = {
          "-i", reads.getPath(),
          "-t", template.getPath(),
          "-o", out.getPath(),
          "-E", "100%", "-m", "0",
          "--sam-rg", header.getPath(),
          "--mask", "cgmaska15b1",
          "--" + MapFlags.SAM_FLAG,
          "--" + MapFlags.DONT_UNIFY_FLAG
      };

      assertEquals(0, map.mainInit(args, TestUtils.getNullOutputStream(), System.err));

      final String mated = FileHelper.gzFileToString(new File(out, "mated.sam.gz"));
      //System.err.println(mated);
      TestUtils.containsAll(mated, "@HD\tVN:1.4\tSO:coordinate",
          "@RG\tID:L23\tPL:COMPLETE\tSM:NA123",
          "@PG\tID:rtg",
          "@CO\tTEMPLATE-SDF-ID:",
          "@CO\tREAD-SDF-ID:",
          "@SQ\tSN:template\tLN:160");
      final String unmated = FileHelper.gzFileToString(new File(out, "unmated.sam.gz"));
      //      System.err.println(unmated);
      TestUtils.containsAll(unmated, "@HD\tVN:1.4\tSO:coordinate",
          "@RG\tID:L23\tPL:COMPLETE\tSM:NA123",
          "@PG\tID:rtg",
          "@CO\tTEMPLATE-SDF-ID:",
          "@CO\tREAD-SDF-ID:",
          "@SQ\tSN:template\tLN:160",
          //old cg produces:            "AS:i:4\tNM:i:3\tGS:Z:TC\tGC:Z:4S1G29S\tGQ:Z:#\tRG:Z:L23"
          "AS:i:2\tXU:Z:5T20=4N2=2X6=\tXR:Z:TACGTCA\tRG:Z:L23\tIH:i:1\tNH:i:1"
          );
      final String usageLog = map.usageLog();
      //System.err.println(usageLog);
      TestUtils.containsAll(usageLog, "[Usage beginning module=cgmap runId=", ", Usage end module=cgmap runId=", " metric=" + LENGTH + " success=true]");
    } finally {
      assertTrue(FileHelper.deleteAll(outer));
    }
  }


  public void testSamRGErrors() throws IOException {
    Diagnostic.setLogStream();
    final File outer = FileUtils.createTempDir("cgmap", "end2end");
    try {
      final File template = new File(outer, "template");
      final File reads = new File(outer, "Reads");
      final File left = new File(reads, "left");
      final File right = new File(reads, "right");
      final File out = new File(outer, "out");

      ReaderTestUtils.getReaderDNA(TEMPLATE, template, null).close();
      ReaderTestUtils.getReaderDNAFastqCG(READ_LEFT, left, PrereadArm.LEFT).close();
      ReaderTestUtils.getReaderDNAFastqCG(READ_RIGHT, right, PrereadArm.RIGHT).close();

      final CgMapCli map = new CgMapCli();

      final String[] args = {
          "-i", reads.getPath(),
          "-t", template.getPath(),
          "-o", out.getPath(),
          "-E", "100%", "-m", "0",
          "--sam-rg", "foo"
      };
      final MemoryPrintStream ps = new MemoryPrintStream();
      assertEquals(1, map.mainInit(args, TestUtils.getNullOutputStream(), ps.printStream()));
      assertTrue(ps.toString(), ps.toString().contains("No read group information present in the input string"));
      final String usageLog = map.usageLog();
      //System.err.println(usageLog);
      TestUtils.containsAll(usageLog, "[Usage beginning module=cgmap runId=", ", Usage end module=cgmap runId=", " metric=0 success=false]");

      final File header = new File(outer, "header");
      FileUtils.stringToFile("", header);
      final String[] args2 = {
          "-i", reads.getPath(),
          "-t", template.getPath(),
          "-o", out.getPath(),
          "-E", "100%", "-m", "0",
          "--sam-rg", header.getPath()
      };
      final MemoryPrintStream ps2 = new MemoryPrintStream();
      assertEquals(1, map.mainInit(args2, TestUtils.getNullOutputStream(), ps2.printStream()));
      assertTrue(ps2.toString(), ps2.toString().contains("No read group information present in the input file \"" + header.getPath() + "\", please provide file with single read group line"));


      final File header2 = new File(outer, "header2");
      FileUtils.stringToFile("@RG\tID:L23\tSM:NA123" + "\n" + "@RG\tID:L43\tSM:NA123", header2);
      final String[] args3 = {
          "-i", reads.getPath(),
          "-t", template.getPath(),
          "-o", out.getPath(),
          "-E", "100%", "-m", "0",
          "--sam-rg", header2.getPath()
      };
      final MemoryPrintStream ps3 = new MemoryPrintStream();
      assertEquals(1, map.mainInit(args3, TestUtils.getNullOutputStream(), ps3.printStream()));
      assertTrue(ps3.toString(), ps3.toString().contains("Multiple read group information present in the input file \"" + header2.getPath() + "\", please provide file with single read group line"));

    } finally {
      assertTrue(FileHelper.deleteAll(outer));
    }
  }


  private static final String CG_QUALITY = "+" + LS + "00000" + "0000000000" + "0000000000" + "0000000000" + LS;
  public void testCGMapbug1144() throws IOException {
    //set showBug to true to get the failure for the bug 1144
    final boolean showBug = false;
    Diagnostic.setLogStream();
    final String templateSeq = ">template" + LS + "GACCATCAGGACAAACACATGGATACATGGAGGGGAACACACACAC" + LS;
    //                                               CAATCAGGAC      GTGGATACATGGAGGGGAAC
    //                                                                                 ACA AC
    //System.err.println(templateSeq);
    final PrintStream err = System.err;
    final File parent = FileUtils.createTempDir("testbug1144", "overlap2Same");
    try {
      if (!showBug) {
        System.setErr(TestUtils.getNullPrintStream());
      }
      final File reads = new File(parent, "reads");
      final File left = new File(reads, "left");
      final File right = new File(reads, "right");
      final String left0 =  "AAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA"; //Random seq
      ReaderTestUtils.getReaderDNAFastqCG("@reL" + LS + left0 + CG_QUALITY, left, PrereadArm.LEFT);
      final String right0 = "CAATCAGGACGTGGATACATGGAGGGGAACACAAC"; //Buggy seq
      ReaderTestUtils.getReaderDNAFastqCG("@reR" + LS + right0 + CG_QUALITY, right, PrereadArm.RIGHT);

      final File template = new File(parent, "template");
      ReaderTestUtils.getReaderDNA(templateSeq, template, null);
      final File output = new File(parent, "out");
      final CgMapCli cgmap = new CgMapCli();
      final MemoryPrintStream mps = new MemoryPrintStream();
      final int code = cgmap.mainInit(new String[] {"-o", output.getPath(), "-i", reads.getPath(), "-t", template.getPath(), "-E", "10", "-Z", "-T", "1", "--sam", "--no-merge"}, mps.outputStream(), mps.printStream());
      assertEquals(mps.toString(), 0, code);
      if (showBug) {
        final String str = dumpSam(new File(output, "unmated.sam"));
        //System.err.println(str);
        //not clear what the real alignment should be - the "good" result below assumes an overlap of 2
        assertTrue(str.contains("GC:Z:28S2G3S")); //bug has GC:Z:30S1G3S
      }
    } finally {
      System.setErr(err);
      //System.err.println(parent.getPath());
      FileHelper.deleteAll(parent);
    }
  }

  String dumpSam(final File samFile) throws IOException {
    final StringBuilder sb = new StringBuilder();
    final RecordIterator<SAMRecord> iter = new SkipInvalidRecordsIterator(samFile);
    while (iter.hasNext()) {
      final SAMRecord sam = iter.next();
      sb.append(sam.getSAMString().trim()).append(LS);
    }
    iter.close();
    return sb.toString();
  }
}
