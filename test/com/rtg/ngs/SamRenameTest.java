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


import java.io.File;
import java.io.IOException;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SdfId;
import com.rtg.reader.SequencesReader;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.diagnostic.CliDiagnosticListener;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.DiagnosticListener;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.diagnostic.SlimException;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

/**
 */
public class SamRenameTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new SamRename();
  }

  /**
   * @param args command line arguments
   */
  public static void main(final String[] args) {
    junit.textui.TestRunner.run(SamRenameTest.class);
  }

  /**
   * Test of moduleName method, of class CgReadSimulator.
   */
  public void testModuleName() {
    final String expResult = "samrename";
    final String result = new SamRename().moduleName();
    assertEquals(expResult, result);
  }

  /**
   * Test of mainInit method, of class SnpStats.
   * @throws IOException if an IO error occurs
   */
  public void testFlags() throws IOException {
    final File tempDir = FileHelper.createTempDirectory();
    try {
      final File samfile = new File(tempDir, "somesamfile");
      //      final File renamedfile = new File(tempDir, "renamedfile");
      //      final File readsdir = new File(tempDir, "detected");
      checkMainInitBadFlags();
      checkMainInitBadFlags(samfile.getAbsolutePath());
    } finally {
      assertTrue(FileHelper.deleteAll(tempDir));
    }
  }

  /**
   * Test of usage message
   */
  public void testUsage() {
    checkHelp("SDF for the reads in the SAM file",
              "input SAM file",
              "print help on command-line flag usage",
              "renamed output SAM file");
  }

  private static final String LS = StringUtils.LS; // Interesting, normally we just use \n for SAM, but we should test readability of \r\n too...
  private static final String TB = "\t";
  private static final String SAM = ""
    + "@HD" + TB + "VN:1.4" + TB + "SO:coordinate" + LS
    + "@SQ" + TB + "SN:g1" + TB + "LN:20" + LS
    + "0" + TB + "0" + TB + "g1" + TB +  "3" + TB + "255" + TB + "8M" + TB + "*" + TB + "0" + TB + "0" + TB + "ATCGACTG" + TB + "*" + TB + "AS:i:0" + LS
    + "0" + TB + "0" + TB + "g1" + TB +  "3" + TB + "255" + TB + "8M" + TB + "*" + TB + "0" + TB + "0" + TB + "ATCGACTG" + TB + "*" + TB + "AS:i:0" + LS
    + "0" + TB + "0" + TB + "g1" + TB +  "5" + TB + "255" + TB + "8M" + TB + "*" + TB + "0" + TB + "0" + TB + "CGACTGTT" + TB + "*" + TB + "AS:i:1" + LS
    + "0" + TB + "0" + TB + "g1" + TB +  "6" + TB + "255" + TB + "8M" + TB + "*" + TB + "0" + TB + "0" + TB + "GACTGCTC" + TB + "*" + TB + "AS:i:1" + LS
    + "0" + TB + "0" + TB + "g1" + TB +  "6" + TB + "255" + TB + "8M" + TB + "*" + TB + "0" + TB + "0" + TB + "GACTGCTC" + TB + "*" + TB + "AS:i:1" + LS
    + "0" + TB + "0" + TB + "g1" + TB + "11" + TB + "255" + TB + "8M" + TB + "*" + TB + "0" + TB + "0" + TB + "TTCAGCTA" + TB + "*" + TB + "AS:i:1" + LS
    ;

  /** Test for an error that will be picked up during params object construction. */
  public void checkParamsError(final String[] args0, final String exp) throws InvalidParamsException {
    TestUtils.containsAll(checkHandleFlagsErr(args0), exp);
  }

  // commandline errors
  public void testCliError1() throws InvalidParamsException, IOException {
    checkParamsError(new String[] {"notexists"}, "You must provide a value for -i SDF");
  }

  // commandline errors
  public void testCliError2() throws InvalidParamsException, IOException {
    final File readsDir = FileUtils.createTempDir("samrenametest", "reads");
    try {
      ReaderTestUtils.getReaderDNA(">r1" + LS + "acgt", readsDir, null).close();
      checkParamsError(new String[] {"notexists", "-i", readsDir.getPath()}, "Input file \"notexists\" doesn't exist");
    } finally {
      assertTrue(FileHelper.deleteAll(readsDir));
    }
  }

  // commandline errors
  public void testCliError3() throws InvalidParamsException, IOException {
    final File readsDir = FileUtils.createTempDir("samrenametest", "reads");
    ReaderTestUtils.getReaderDNA(">r1" + LS + "acgt", readsDir, null).close();
    final File tempDir = FileUtils.createTempDir("samrenametest", "checkparamserror");
    try {
      final File map = new File(tempDir, "map");
      FileUtils.stringToFile(SAM, map);
      final File mapRename = new File(tempDir, "map_rename");
      checkMainInitWarn(map.getPath(), "-i", readsDir.getPath(),
                        "-o", mapRename.getPath(), "--no-gzip");
      assertTrue(mapRename.exists());
    } finally {
      assertTrue(FileHelper.deleteAll(readsDir));
      assertTrue(FileHelper.deleteAll(tempDir));
    }
  }

  // commandline errors
  public void testCliError3b() throws InvalidParamsException, IOException {
    final File readsDir = FileUtils.createTempDir("samrenametest", "reads");
    ReaderTestUtils.getReaderDNA(">r1" + LS + "acgt", readsDir, null).close();
    final File tempDir = FileUtils.createTempDir("samrenametest", "checkparamserror");
    try {
      final File map = new File(tempDir, "map");
      FileUtils.stringToFile(SAM, map);
      final File mapRename = new File(tempDir, "map_rename");
      checkMainInitWarn(map.getPath(), "-i", readsDir.getPath(),
                        "--output", mapRename.getPath(), "--no-gzip");
      assertTrue(mapRename.exists());
    } finally {
      assertTrue(FileHelper.deleteAll(readsDir));
      assertTrue(FileHelper.deleteAll(tempDir));
    }
  }

  public void testCliError4() throws InvalidParamsException, IOException {
    final File readsDir = FileUtils.createTempDir("samrenametest", "reads");
    ReaderTestUtils.getReaderDNA(">r1" + LS + "acgt", readsDir, null).close();
    final File tempDir = FileUtils.createTempDir("samrenametest", "checkparamserror");
    try {
      final File map = new File(tempDir, "map.sam");
      FileUtils.stringToFile(SAM, map);
      checkMainInitWarn(map.getPath(), "-i", readsDir.getPath(), "--no-gzip");
      final File mapRename = new File(tempDir, "map_rename.sam");
      assertTrue(mapRename.exists());
      TestUtils.containsAll(FileUtils.fileToString(mapRename), "@PG\tID:rtg");
    } finally {
      assertTrue(FileHelper.deleteAll(readsDir));
      assertTrue(FileHelper.deleteAll(tempDir));
    }
  }

  public void testCliError4Gzipped() throws InvalidParamsException, IOException {
    final File readsDir = FileUtils.createTempDir("samrenametest", "reads");
    ReaderTestUtils.getReaderDNA(">r1" + LS + "acgt", readsDir, null).close();
    final File tempDir = FileUtils.createTempDir("samrenametest", "checkparamserror");
    try {
      final File map = new File(tempDir, "map.sam.gz");
      FileHelper.stringToGzFile(SAM, map);

      checkMainInitWarn(map.getPath(), "-i", readsDir.getPath());
      final File mapRename = new File(tempDir, "map_rename.sam.gz");
      assertTrue(mapRename.exists());
    } finally {
      assertTrue(FileHelper.deleteAll(readsDir));
      assertTrue(FileHelper.deleteAll(tempDir));
    }
  }

  private static final String SAM2_HEAD = ""
    + "@HD" + TB + "VN:1.4" + TB + "SO:coordinate" + "\n"
    + "@SQ" + TB + "SN:g1" + TB + "LN:20" + "\n";
  private static final String SAM2_DATA = ""
    + "0" + TB + "0" + TB + "g1" + TB +  "3" + TB + "255" + TB + "8M" + TB + "*" + TB + "0" + TB + "0" + TB + "ATCGACTG" + TB + "*" + TB + "AS:i:0" + "\n"
    + "0" + TB + "0" + TB + "g1" + TB +  "3" + TB + "255" + TB + "8M" + TB + "*" + TB + "0" + TB + "0" + TB + "ATCGACTG" + TB + "*" + TB + "AS:i:0" + "\n"
    + "0" + TB + "0" + TB + "g1" + TB +  "5" + TB + "255" + TB + "8M" + TB + "*" + TB + "0" + TB + "0" + TB + "CGACTGTT" + TB + "*" + TB + "AS:i:1" + "\n"
    + "0" + TB + "0" + TB + "g1" + TB +  "6" + TB + "255" + TB + "8M" + TB + "*" + TB + "0" + TB + "0" + TB + "GACTGCTC" + TB + "*" + TB + "AS:i:1" + "\n"
    + "0" + TB + "0" + TB + "g1" + TB +  "6" + TB + "255" + TB + "8M" + TB + "*" + TB + "0" + TB + "0" + TB + "GACTGCTC" + TB + "*" + TB + "AS:i:1" + "\n"
    + "0" + TB + "0" + TB + "g1" + TB + "11" + TB + "255" + TB + "8M" + TB + "*" + TB + "0" + TB + "0" + TB + "TTCAGCTA" + TB + "*" + TB + "AS:i:1" + "\n"
    ;

  private static final String SAM2_DATA_EXP = ""
    + "r1" + TB + "0" + TB + "g1" + TB +  "3" + TB + "255" + TB + "8M" + TB + "*" + TB + "0" + TB + "0" + TB + "ATCGACTG" + TB + "*" + TB + "AS:i:0" + "\n"
    + "r1" + TB + "0" + TB + "g1" + TB +  "3" + TB + "255" + TB + "8M" + TB + "*" + TB + "0" + TB + "0" + TB + "ATCGACTG" + TB + "*" + TB + "AS:i:0" + "\n"
    + "r1" + TB + "0" + TB + "g1" + TB +  "5" + TB + "255" + TB + "8M" + TB + "*" + TB + "0" + TB + "0" + TB + "CGACTGTT" + TB + "*" + TB + "AS:i:1" + "\n"
    + "r1" + TB + "0" + TB + "g1" + TB +  "6" + TB + "255" + TB + "8M" + TB + "*" + TB + "0" + TB + "0" + TB + "GACTGCTC" + TB + "*" + TB + "AS:i:1" + "\n"
    + "r1" + TB + "0" + TB + "g1" + TB +  "6" + TB + "255" + TB + "8M" + TB + "*" + TB + "0" + TB + "0" + TB + "GACTGCTC" + TB + "*" + TB + "AS:i:1" + "\n"
    + "r1" + TB + "0" + TB + "g1" + TB + "11" + TB + "255" + TB + "8M" + TB + "*" + TB + "0" + TB + "0" + TB + "TTCAGCTA" + TB + "*" + TB + "AS:i:1" + "\n"
    ;

  private void checkSdfId(File dir, File readsDir, String expOut, String expErr, boolean succeed, SdfId samSdf) throws IOException {
    CommandLine.clearCommandArgs();
    final String commentHeader = !samSdf.available() ? "" : "@CO\tREAD-SDF-ID:" + samSdf + "\n";

    final MemoryPrintStream out = new MemoryPrintStream();
    final MemoryPrintStream err = new MemoryPrintStream();
    final DiagnosticListener dl = new CliDiagnosticListener(err.printStream(), out.printStream());
    Diagnostic.addListener(dl);
    try {
      final File outSam = new File(dir, "outSam.sam");
      final File inSam = new File(dir, "inSam.sam");
      FileUtils.stringToFile(SAM2_HEAD + commentHeader + SAM2_DATA, inSam);
      final SamRename sr = new SamRename();
      try {
        sr.renameSam(readsDir, inSam, outSam, LongRange.NONE);
        if (!succeed) {
          fail("shouldn't have passed");
        }
        final String outSamStr = FileUtils.fileToString(outSam);
        assertEquals(SAM2_HEAD + commentHeader + SAM2_DATA_EXP, outSamStr.replaceAll("@PG\t[^\\n]*\\n", ""));
        assertTrue(outSamStr.contains("@PG\tID:rtg"));
        assertTrue(outSamStr.contains("CL:Internal"));
      } catch (final SlimException e) {
        e.printErrorNoLog();
        if (succeed) {
          fail("expected to pass");
        }
      } catch (final RuntimeException e) {
        if (succeed) {
          fail("expected to pass");
        }
      }
    } finally {
      Diagnostic.removeListener(dl);
    }
    assertEquals(expOut, out.toString());
    assertEquals(expErr, err.toString());
  }

  public void testSdfIdCheck() throws IOException {
    final File dir = FileUtils.createTempDir("samrename", "sdfidcheck");
    try {
      final File readsDir = new File(dir, "reads");
      final SequencesReader reader = ReaderTestUtils.getReaderDNA(">r1" + LS + "acgt", readsDir, null);
      try {
        final SdfId sdfId = reader.getSdfId();
        checkSdfId(dir, readsDir, "", "", true, sdfId);
      } finally {
        reader.close();
      }
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  public void testSdfIdCheckWrong() throws IOException {
    final File dir = FileUtils.createTempDir("samrename", "sdfidcheck");
    try {
      final File readsDir = new File(dir, "reads");
      final SequencesReader reader = ReaderTestUtils.getReaderDNA(">r1" + LS + "acgt", readsDir, null);
      try {
        checkSdfId(dir, readsDir, "", "Error: SDF-ID of given SDF does not match SDF used during mapping." + LS, false, new SdfId());
      } finally {
        reader.close();
      }
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  public void testSdfIdCheckUnknown() throws IOException {
    final File dir = FileUtils.createTempDir("samrename", "sdfidcheck");
    try {
      final File readsDir = new File(dir, "reads");
      final SequencesReader reader = ReaderTestUtils.getReaderDNA(">r1" + LS + "acgt", readsDir, null);
      try {
        checkSdfId(dir, readsDir, "", "No READ-SDF-ID found in SAM header, unable to verify read-id correctness." + LS, true, new SdfId(0));
      } finally {
        reader.close();
      }
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  private static final String READS_10 = "\n"
          + ">zero\nacgttgca"
          + ">one\nacgttgca"
          + ">two\nacgttgca"
          + ">three\nacgttgca"
          + ">four\nacgttgca"
          + ">five\nacgttgca"
          + ">six\nacgttgca"
          + ">seven\nacgttgca"
          + ">eight\nacgttgca"
          + ">nine\nacgttgca";
  private static final String SAM3_DATA = ""
    + "2" + TB + "0" + TB + "g1" + TB +  "3" + TB + "255" + TB + "8M" + TB + "*" + TB + "0" + TB + "0" + TB + "ATCGACTG" + TB + "*" + TB + "AS:i:0" + "\n"
    + "3" + TB + "0" + TB + "g1" + TB +  "3" + TB + "255" + TB + "8M" + TB + "*" + TB + "0" + TB + "0" + TB + "ATCGACTG" + TB + "*" + TB + "AS:i:0" + "\n"
    + "4" + TB + "0" + TB + "g1" + TB +  "5" + TB + "255" + TB + "8M" + TB + "*" + TB + "0" + TB + "0" + TB + "CGACTGTT" + TB + "*" + TB + "AS:i:1" + "\n"
    + "5" + TB + "0" + TB + "g1" + TB +  "6" + TB + "255" + TB + "8M" + TB + "*" + TB + "0" + TB + "0" + TB + "GACTGCTC" + TB + "*" + TB + "AS:i:1" + "\n"
    + "6" + TB + "0" + TB + "g1" + TB +  "6" + TB + "255" + TB + "8M" + TB + "*" + TB + "0" + TB + "0" + TB + "GACTGCTC" + TB + "*" + TB + "AS:i:1" + "\n"
    + "7" + TB + "0" + TB + "g1" + TB + "11" + TB + "255" + TB + "8M" + TB + "*" + TB + "0" + TB + "0" + TB + "TTCAGCTA" + TB + "*" + TB + "AS:i:1" + "\n"
    + "8" + TB + "0" + TB + "g1" + TB + "11" + TB + "255" + TB + "8M" + TB + "*" + TB + "0" + TB + "0" + TB + "TTCAGCTA" + TB + "*" + TB + "AS:i:1" + "\n"
    ;

  private static final String SAM3_DATA_EXP = ""
    + "two" + TB + "0" + TB + "g1" + TB +  "3" + TB + "255" + TB + "8M" + TB + "*" + TB + "0" + TB + "0" + TB + "ATCGACTG" + TB + "*" + TB + "AS:i:0" + "\n"
    + "three" + TB + "0" + TB + "g1" + TB +  "3" + TB + "255" + TB + "8M" + TB + "*" + TB + "0" + TB + "0" + TB + "ATCGACTG" + TB + "*" + TB + "AS:i:0" + "\n"
    + "four" + TB + "0" + TB + "g1" + TB +  "5" + TB + "255" + TB + "8M" + TB + "*" + TB + "0" + TB + "0" + TB + "CGACTGTT" + TB + "*" + TB + "AS:i:1" + "\n"
    + "five" + TB + "0" + TB + "g1" + TB +  "6" + TB + "255" + TB + "8M" + TB + "*" + TB + "0" + TB + "0" + TB + "GACTGCTC" + TB + "*" + TB + "AS:i:1" + "\n"
    + "six" + TB + "0" + TB + "g1" + TB +  "6" + TB + "255" + TB + "8M" + TB + "*" + TB + "0" + TB + "0" + TB + "GACTGCTC" + TB + "*" + TB + "AS:i:1" + "\n"
    + "seven" + TB + "0" + TB + "g1" + TB + "11" + TB + "255" + TB + "8M" + TB + "*" + TB + "0" + TB + "0" + TB + "TTCAGCTA" + TB + "*" + TB + "AS:i:1" + "\n"
    + "eight" + TB + "0" + TB + "g1" + TB + "11" + TB + "255" + TB + "8M" + TB + "*" + TB + "0" + TB + "0" + TB + "TTCAGCTA" + TB + "*" + TB + "AS:i:1" + "\n"
    ;

  public void testMoreThanOne() throws IOException {
    try (TestDirectory dir = new TestDirectory()) {
      final File readsDir = new File(dir, "reads");
      final File samIn = FileUtils.stringToFile(SAM2_HEAD + SAM3_DATA, new File(dir, "samIn.sam"));
      final File samOut = new File(dir, "samOut.sam");
      ReaderTestUtils.getReaderDNA(READS_10, readsDir, null).close();
      final SamRename sr = new SamRename();
      sr.renameSam(readsDir, samIn, samOut, new LongRange(2, 9));
      assertEquals(SAM2_HEAD + SAM3_DATA_EXP, FileUtils.fileToString(samOut).replaceAll("@PG[^\\n]*\\n", ""));
      assertTrue(samOut.delete());
      sr.renameSam(readsDir, samIn, samOut, LongRange.NONE);
      assertEquals(SAM2_HEAD + SAM3_DATA_EXP, FileUtils.fileToString(samOut).replaceAll("@PG[^\\n]*\\n", ""));
      assertTrue(samOut.delete());
      sr.renameSam(readsDir, samIn, samOut, new LongRange(0, 50));
      assertEquals(SAM2_HEAD + SAM3_DATA_EXP, FileUtils.fileToString(samOut).replaceAll("@PG[^\\n]*\\n", ""));
    }
  }

  public void testErrors() throws IOException {
    try (TestDirectory dir = new TestDirectory()) {
      final File readsDir = new File(dir, "reads");
      final File samIn = FileUtils.stringToFile(SAM2_HEAD + SAM3_DATA, new File(dir, "samIn.sam"));
      final File samOut = new File(dir, "samOut.sam");
      ReaderTestUtils.getReaderDNA(READS_10, readsDir, null).close();
      final SamRename sr = new SamRename();
      try {
        sr.renameSam(readsDir, samIn, samOut, new LongRange(3, 9));
        fail();
      } catch (NoTalkbackSlimException e) {
        assertTrue(e.getMessage().contains("lower"));
      }
      try {
        sr.renameSam(readsDir, samIn, samOut, new LongRange(2, 8));
        fail();
      } catch (NoTalkbackSlimException e) {
        assertTrue(e.getMessage().contains("higher"));
      }
      try {
        sr.renameSam(readsDir, samIn, samOut, new LongRange(3, 6));
        fail();
      } catch (NoTalkbackSlimException e) {
        assertTrue(e.getMessage().contains("lower"));
      }
    }
  }
}
