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
package com.rtg.launcher;


import static com.rtg.launcher.BuildCommon.RESOURCE;
import static com.rtg.launcher.CommonFlags.INT;
import static com.rtg.util.cli.CommonFlagCategories.REPORTING;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.util.List;

import com.rtg.calibrate.Recalibrate;
import com.rtg.ngs.OutputFilter;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.IntegerOrPercentage;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.Flag;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;
import com.rtg.util.test.MockEventListener;

import junit.framework.TestCase;


/**
 */
public class CommonFlagsTest extends TestCase {

  private static final String LONG_READ_1 = "CATAATGACGGCTGGGCTACTGGACATATGTACGCGGTCTCCGGGAACAAGCAGCATGGAGTTTCCCCTGACGTATCGGTGATGTGATTGACATAACGATCAGATTTCAAAAGGAGTTCGCGCATTCCAGAGGACGCTATGCACGTTGGT";
  private static final String LONG_READ_2 = "ACAAGTACAAGGAGGCGGACTATCGCACTCCAAGTTAGCCGGTTTGAACATAGAAGCTCTCCCGAGCTCGCGCTATCCATCGATGTGGATGAGGAGTCAGCAATAAACTCGCTGATAGCCTGATCCAGTCAACTTTATCCATGTCAGGCA";

  private static final String TEMPLATE_LONG = LONG_READ_1 + "atcgaggtcatctagcagcatcatcgacttatcgacatctacgatcgagcgactatcgactatg" + LONG_READ_2;

  public void testInitMaskFlag() {
    final CFlags flags = new CFlags("blah", TestUtils.getNullPrintStream(), TestUtils.getNullPrintStream());
    CommonFlags.initMaskFlags(flags);

    assertNotNull(flags.getFlag("word"));
    assertNotNull(flags.getFlag("substitutions"));
    assertNotNull(flags.getFlag("indels"));
    assertNotNull(flags.getFlag("indel-length"));

    TestUtils.containsAll(flags.getUsageString(),
        " -w,", "word size",
        " -a,", "guaranteed minimum number of substitutions",
        " -b,", "guaranteed minimum number of indels",
        " -c,", "guaranteed number of positions"
    );
  }

  public void testSharedFlag() {
    final CFlags flags = new CFlags("blah", TestUtils.getNullPrintStream(), TestUtils.getNullPrintStream());
    CommonFlags.initSharedFlags(flags);

    assertNotNull(flags.getFlag(CommonFlags.STEP_FLAG));
    assertNotNull(flags.getFlag(CommonFlags.THREADS_FLAG));
    assertNotNull(flags.getFlag(CommonFlags.NO_GZIP));

    TestUtils.containsAll(flags.getUsageString(500),
        "-s", "step size (Default is word size)",
        RESOURCE.getString("REPEAT_FREQUENCY_DESC"),
        "-T", "number of threads. Defaults to the number of available cores",
        "-Z", "do not gzip the output"
        );
  }

  public void testValidateMask() {
    valMaskInvalidValue(0, 0, 0, 0, 0, "w", 0, 1, true);
    valMaskInvalidValue(0, 1, -2, 0, 0, "a", -2, 0, true);
    valMaskInvalidValue(0, 1, 2, -3, 0, "b", -3, 0, true);
    valMaskInvalidValue(0, 1, 2, 3, 0, "c", 0, 1, true);

    final ByteArrayOutputStream baos = new ByteArrayOutputStream();
    final PrintStream ps = new PrintStream(baos);
    Diagnostic.setLogStream(ps);
    try {
      CommonFlags.validateMaskParams(5, 6, 0, 0, 1);
      fail();
    } catch (final InvalidParamsException ipe) {
      assertEquals(ErrorType.WORD_NOT_LESS_READ, ipe.getErrorType());
//      ps.flush();
//      assertTrue(baos.toString().contains("The word length \"6\" should be less than the read length \"5\""));
    }
    Diagnostic.setLogStream();

    valMaskInvalidValue(5, 1, 2, 3, 6, "c", 6, 4, false);
    valMaskInvalidValue(5, 1, 2, 3, 0, "c", 0, 1, true);

    valMaskInvalidValue(5, 1, 2, 6, 3, "b", 6, 4, false);
    valMaskInvalidValue(5, 1, 6, 3, 3, "a", 6, 4, false);
  }

  private void valMaskInvalidValue(int readlen, int w, int subs, int indels, int indellen, String badflag, int badvalue, int expected, boolean max) {
    final ByteArrayOutputStream baos = new ByteArrayOutputStream();
    final PrintStream ps = new PrintStream(baos);
    Diagnostic.setLogStream(ps);
    try {
      CommonFlags.validateMaskParams(readlen, w, subs, indels, indellen);
      fail();
    } catch (final InvalidParamsException ipe) {
      assertEquals(max ? ErrorType.INVALID_MIN_INTEGER_FLAG_VALUE : ErrorType.INVALID_MAX_INTEGER_FLAG_VALUE, ipe.getErrorType());
      ps.flush();
//      final String out = baos.toString();
//      final String exp = "The specified flag \"-" + badflag + "\" has invalid value \"" + badvalue + "\". It should be "
//              + (max ? "greater than or equal to " : "less than or equal to ")
//              + "\"" + expected + "\"";
//
//      assertTrue(ipe.getMessage().contains(exp));
    }
    Diagnostic.setLogStream();
  }

  public void testOutputFilter() {
    final CFlags flags = new CFlags();

    flags.registerOptional('f', CommonFlags.OUTPUT_FILTER_FLAG, String.class, "name", "output filter", "none");
    flags.registerOptional(CommonFlags.TOPN_RESULTS_FLAG, Integer.class, INT, "set the number of results per read for topn. Allowed values are between 1 and 255", 5).setCategory(REPORTING);

    assertEquals(OutputFilter.NONE, CommonFlags.filter(flags));

    assertTrue(flags.setFlags("-f", "PROTEIN_TOPN", "--" + CommonFlags.TOPN_RESULTS_FLAG, "0"));
    assertEquals(OutputFilter.PROTEIN_TOPN, CommonFlags.filter(flags));
  }

  public void testNegativeRepeatFreq() {
    final MockEventListener ev = new MockEventListener();
    Diagnostic.addListener(ev);
    final CFlags flags = new CFlags();
    flags.registerOptional('r', CommonFlags.REPEAT_FREQUENCY_FLAG, Integer.class, "int", RESOURCE.getString("REPEAT_FREQUENCY_DESC"), 1000);
    try {
      flags.setFlags("-r", "-1");
      assertFalse(CommonFlags.checkRepeatFrequency(flags));
      assertTrue(ev.compareErrorMessage("Error: The specified flag \"--" + CommonFlags.REPEAT_FREQUENCY_FLAG + "\" has invalid value \"" + -1 + "\". It should be greater than or equal to \"1\"."));

      flags.setFlags("-r", "0");
      assertFalse(CommonFlags.checkRepeatFrequency(flags));
      assertTrue(ev.compareErrorMessage("Error: The specified flag \"--" + CommonFlags.REPEAT_FREQUENCY_FLAG + "\" has invalid value \"" + 0 + "\". It should be greater than or equal to \"1\"."));
    } finally {
      Diagnostic.removeListener(ev);
    }
  }


  public void testPercentRepeatFreq() {
    final MockEventListener ev = new MockEventListener();
    Diagnostic.addListener(ev);
    final CFlags flags = new CFlags();
    flags.registerOptional('r', CommonFlags.REPEAT_FREQUENCY_FLAG, IntegerOrPercentage.class, "int", RESOURCE.getString("REPEAT_FREQUENCY_DESC"), new IntegerOrPercentage(20));
    try {
      flags.setFlags("-r", "-1");
      assertFalse(CommonFlags.checkPercentRepeatFrequency(flags));
      assertTrue(ev.compareErrorMessage("Error: The specified flag \"--" + CommonFlags.REPEAT_FREQUENCY_FLAG + "\" has invalid value \"" + -1 + "\". It should be greater than or equal to \"1\"."));

      flags.setFlags("-r", "0");
      assertFalse(CommonFlags.checkPercentRepeatFrequency(flags));
      assertTrue(ev.compareErrorMessage("Error: The specified flag \"--" + CommonFlags.REPEAT_FREQUENCY_FLAG + "\" has invalid value \"" + 0 + "\". It should be greater than or equal to \"1\"."));

      flags.setFlags("-r", "100001");
      assertFalse(CommonFlags.checkPercentRepeatFrequency(flags));
      assertTrue(ev.compareErrorMessage("Error: The specified flag \"--" + CommonFlags.REPEAT_FREQUENCY_FLAG + "\" has invalid value \"" + 100001 + "\". It should be less than or equal to \"100000\"."));

      flags.setFlags("-r", "-1%");
      assertFalse(CommonFlags.checkPercentRepeatFrequency(flags));
      assertTrue(ev.compareErrorMessage("Error: The specified flag \"--" + CommonFlags.REPEAT_FREQUENCY_FLAG + "\" has invalid value \"" + -1 + "%\". It should be greater than or equal to \"0%\"."));

      flags.setFlags("-r", "101%");
      assertFalse(CommonFlags.checkPercentRepeatFrequency(flags));
      assertTrue(ev.compareErrorMessage("Error: The specified flag \"--" + CommonFlags.REPEAT_FREQUENCY_FLAG + "\" has invalid value \"" + 101 + "%\". It should be less than or equal to \"100%\"."));

      flags.setFlags("-r", "1");
      assertTrue(CommonFlags.checkPercentRepeatFrequency(flags));
    } finally {
      Diagnostic.removeListener(ev);
    }
  }

  public void testMapIOFlags() {
    final CFlags flags = new CFlags("blah", TestUtils.getNullPrintStream(), TestUtils.getNullPrintStream());
    CommonFlags.initMapIOFlags(flags);
    assertNotNull(flags.getFlag("input"));
    assertNotNull(flags.getFlag("output"));
    assertNotNull(flags.getFlag("template"));

    TestUtils.containsAll(flags.getUsageString(),
        " -i,", "SDF containing reads to map",
        " -o,", "directory for output",
        " -t,", "SDF containing template to map against"
        );
  }

  public void testValidateTemplate() throws IOException {
    final File tempDir = FileUtils.createTempDir("junit", "test");

    final MockEventListener ev = new MockEventListener();
    Diagnostic.addListener(ev);
    final CFlags flags = new CFlags("blah", TestUtils.getNullPrintStream(), TestUtils.getNullPrintStream());
    CommonFlags.initMapIOFlags(flags);
    try {
      flags.setFlags("-t", "template");
      assertFalse(CommonFlags.validateSDF(flags, "template"));
      assertTrue(ev.compareErrorMessage("Error: The specified SDF, \"template\", does not exist."));
      assertTrue(ev.compareErrorMessage("Error: The specified SDF, \"template\", does not exist."));
      final File temp = new File(tempDir, "temp.sdf");
      assertTrue(temp.createNewFile());
      flags.setFlags("-t", temp.getAbsolutePath());
      assertFalse(CommonFlags.validateSDF(flags, "template"));
      assertTrue(ev.compareErrorMessage("Error: The specified file, \"" + temp.getAbsolutePath() + "\", is not an SDF."));
      assertTrue(temp.delete());
      ReaderTestUtils.getReaderDNA(">template\n" + TEMPLATE_LONG, temp, null);
      assertTrue(CommonFlags.validateTemplate(flags));
      assertTrue(FileHelper.deleteAll(temp));

      assertTrue(CommonFlags.validateOutputDirectory(new File(tempDir, "blsdkjrlks")));
      assertFalse(CommonFlags.validateOutputDirectory(tempDir));
    } finally {
      Diagnostic.removeListener(ev);
      assertTrue(FileHelper.deleteAll(tempDir));
    }
  }

  public void testWordSize() {
    final CFlags flags = new CFlags();
    CommonFlags.initWordSize(flags, "");
    assertEquals(5, CommonFlags.getWordSize(flags, 10, 7));
    assertEquals(7, CommonFlags.getWordSize(flags, 15, 7));
    flags.setFlags("--" + CommonFlags.WORDSIZE_FLAG, "8");
    assertEquals(8, CommonFlags.getWordSize(flags, 15, 7));
  }

  public void testThreads() {
    Diagnostic.setLogStream(TestUtils.getNullPrintStream());
    try {
      final CFlags flags = new CFlags("blah", TestUtils.getNullPrintStream(), TestUtils.getNullPrintStream());
      CommonFlags.initSharedFlagsOnly(flags);
      flags.setFlags("-T", "0");
      assertFalse(CommonFlags.validateThreads(flags));
      flags.setFlags("-T", "-1");
      assertFalse(CommonFlags.validateThreads(flags));
      flags.setFlags("-T", Integer.toString(Integer.MAX_VALUE));
      assertFalse(CommonFlags.validateThreads(flags));
      flags.setFlags("-T", "20");
      assertTrue(CommonFlags.validateThreads(flags));
    } finally {
      Diagnostic.setLogStream();
    }
  }

  public void testFileLists() throws IOException {
    final ByteArrayOutputStream baos = new ByteArrayOutputStream();
    Diagnostic.setLogStream(new PrintStream(baos));
    try {
      final String listFlag = "file-list";
      final File tmp = FileHelper.createTempDirectory();
      try {
        final File[] files = new File[10];
        final File listfile = new File(tmp, "file-list");
        try (FileWriter fw = new FileWriter(listfile)) {
          fw.append("# some kind of header to be ignored").append(StringUtils.LS).append("   ").append(StringUtils.LS); // Test comment skipping
          for (int i = 0; i < files.length; i++) {
            files[i] = new File(tmp, "file" + i);
            fw.append(" ").append(files[i].getPath()).append(" ").append(StringUtils.LS); // Test line trimming
          }
        }
        final CFlags flags = new CFlags("blah", TestUtils.getNullPrintStream(), TestUtils.getNullPrintStream());
        flags.registerOptional(listFlag, File.class, "FILE", "files");
        final Flag reads = flags.registerRequired(File.class, "File", "input sam files");
        reads.setMinCount(0);
        reads.setMaxCount(4000);
        final String[] args = new String[files.length];
        for (int i = 0; i < files.length; i++) {
          args[i] = files[i].getPath();
        }
        flags.setFlags(args);
        final CFlags flags2 = new CFlags("blah", TestUtils.getNullPrintStream(), TestUtils.getNullPrintStream());
        flags2.registerOptional(listFlag, File.class, "FILE", "files");
        //final Flag reads2 =
        flags2.registerRequired(File.class, "File", "input sam files");
        reads.setMinCount(0);
        reads.setMaxCount(4000);
        flags2.setFlags("--" + listFlag, listfile.getAbsolutePath());
        try {
          CommonFlags.getFileList(flags, listFlag, null, false);
          fail();
        } catch (final NoTalkbackSlimException e) {
          assertTrue(e.getMessage().contains("There were 10 invalid input file paths"));
        }
        for (int i = 0; i < files.length; i++) {
          assertEquals(i < 5, baos.toString().contains("File not found: \"" + files[i].getPath() + "\""));
        }
        baos.reset();

        for (final File file : files) {
          assertTrue(file.createNewFile());
        }
        CommonFlags.getFileList(flags, listFlag, null, false);
        CommonFlags.getFileList(flags2, listFlag, null, false);
        try {
          CommonFlags.getFileList(flags, listFlag, null, true);
          fail();
        } catch (final NoTalkbackSlimException e) {
          assertTrue(e.getMessage().contains("There were 10 invalid input file paths"));
        }

        for (int i = 0; i < files.length; i++) {

          assertEquals("files[" + i + "]=" + files[i].getPath() + " was " + (i < 5 ? "not " : "") + "contained in the string", i < 5, baos.toString().contains(files[i].getPath() + " is not a valid SDF"));
        }

        for (int i = 0; i < files.length; i++) {
          assertTrue(files[i].delete());
          assertTrue(files[i].mkdir());
        }
        baos.reset();
        try {
          CommonFlags.getFileList(flags, listFlag, null, false);
                  System.err.println(baos.toString());
          fail();
        } catch (final NoTalkbackSlimException e) {
          assertTrue(e.getMessage().contains("There were 10 invalid input file paths"));
        }
        for (int i = 0; i < files.length; i++) {
          assertEquals(i < 5, baos.toString().contains(files[i].getPath() + "\" is not a file"));
        }

        baos.reset();
        try {
          CommonFlags.getFileList(flags2, listFlag, null, false);
          fail();
        } catch (final NoTalkbackSlimException e) {
          assertTrue(e.getMessage(), e.getMessage().contains("There were 10 invalid input file paths"));
        }
        for (int i = 0; i < files.length; i++) {
          assertEquals(i < 5, baos.toString().contains(files[i].getPath() + "\" is not a file"));
        }


        final File f1 = new File(tmp, "f1");
        assertTrue(f1.createNewFile());
        final File f2 = new File(tmp, "f2");
        assertTrue(f2.createNewFile());
        flags.setFlags(f1.getPath(), f2.getPath());
        final List<File> filesout = CommonFlags.getFileList(flags, listFlag, null, false);
        assertNotNull(filesout);
        assertEquals(2, filesout.size());
      } finally {
        assertTrue(FileHelper.deleteAll(tmp));
      }
    } finally {
      Diagnostic.setLogStream();
    }
  }

  public void testCheckFile() {
    final CFlags flags = new CFlags("blah", TestUtils.getNullPrintStream(), TestUtils.getNullPrintStream());
    final Flag inFlag = flags.registerRequired(File.class, "FILE", "");
    inFlag.setMinCount(0);
    inFlag.setMaxCount(Integer.MAX_VALUE);
    flags.registerOptional("input", File.class, "FILE", "i");

    final String[] files = {"bah", "buh"};
    flags.setFlags(files);

    assertFalse(CommonFlags.checkFileList(flags, "blah", null, 1));
    assertEquals("More than 1 input files specified.", flags.getParseMessage());

    flags.setFlags();

    CommonFlags.checkFileList(flags, "blah", "input", 5);
    assertEquals("No input files specified in --blah or --input.", flags.getParseMessage());

    flags.setFlags("abcd", "bcde" + Recalibrate.EXTENSION);
    assertFalse(CommonFlags.checkFileList(flags, "blah", null, 1, false));
    assertTrue(CommonFlags.checkFileList(flags, "blah", null, 1, true));
  }

  public void testReaderRestriction() {
    final CFlags flags = new CFlags("blah", TestUtils.getNullPrintStream(), TestUtils.getNullPrintStream());
    CommonFlags.initReadRange(flags);

    TestUtils.containsAll(flags.getUsageString(),
        "--end-read=INT", "exclusive upper bound on read id",
        "--start-read=INT", "inclusive lower bound on read id"
        );

    LongRange r = CommonFlags.getReaderRestriction(flags);

    assertEquals(-1, r.getStart());
    assertEquals(-1, r.getEnd());

    flags.setFlags("--" + CommonFlags.START_READ_ID, "3", "--" + CommonFlags.END_READ_ID, "5");

    r = CommonFlags.getReaderRestriction(flags);

    assertEquals(3, r.getStart());
    assertEquals(5, r.getEnd());
  }

  public void testValidateSDF() throws Exception {
    final File tmpFile = FileUtils.createTempDir("commonflags", "tmp");
    final MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());
    try {

      final CFlags flags = new CFlags("blah", TestUtils.getNullPrintStream(), TestUtils.getNullPrintStream());
      flags.registerRequired('i', CommonFlags.READS_FLAG, File.class, "SDF", "");
      final File sdf = new File(tmpFile, "sdf");
      flags.setFlags("-i", sdf.getPath());

      assertFalse(CommonFlags.validateReads(flags, true));
      assertTrue(mps.toString(), mps.toString().contains("The specified SDF, \"" + sdf.getPath() + "\", does not exist"));
      mps.reset();

      assertTrue(sdf.createNewFile());

      assertFalse(CommonFlags.validateReads(flags, true));
      assertTrue(mps.toString(), mps.toString().contains("The specified file, \"" + sdf.getPath() + "\", is not an SDF"));
      mps.reset();

      assertTrue(sdf.delete());
      assertTrue(sdf.mkdir());

      assertTrue(mps.toString(), CommonFlags.validateReads(flags, true));
    } finally {
      Diagnostic.setLogStream();
      FileHelper.deleteAll(tmpFile);
    }
  }

  public void testValidateStartEnd() {
    final MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());
    try {

      final CFlags flags = new CFlags("blah", TestUtils.getNullPrintStream(), mps.printStream());
      flags.registerOptional("start", Long.class, "LONG", "s");
      flags.registerOptional("end", Long.class, "LONG", "e");

      flags.setFlags("--start", "-1");
      assertFalse(CommonFlags.validateStartEnd(flags, "start", "end"));
      assertTrue(mps.toString().contains("--start should be positive"));

      mps.reset();
      flags.setFlags("--end", "-1");
      assertFalse(CommonFlags.validateStartEnd(flags, "start", "end"));
      assertTrue(mps.toString().contains("--end should be greater than 0"));

      mps.reset();
      flags.setFlags("--start", "3", "--end", "1");
      assertFalse(CommonFlags.validateStartEnd(flags, "start", "end"));
      assertTrue(mps.toString().contains("--start should be less than --end"));

      mps.reset();
      flags.setFlags("--start", "0", "--end", "" + Long.MAX_VALUE);
      assertFalse(CommonFlags.validateStartEnd(flags, "start", "end"));
      assertTrue(mps.toString().contains("You have specified too many reads, please specify a range of less than " + Integer.MAX_VALUE + " reads."));
    } finally {
      Diagnostic.setLogStream();
    }
  }

  public void testNoSuchAvrModel() {
    final MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());
    try {
      final CFlags flags = new CFlags("blah", TestUtils.getNullPrintStream(), mps.printStream());
      flags.registerOptional(CommonFlags.AVR_MODEL_FILE_FLAG, File.class, "File", "AVR");
      flags.setFlags("--avr-model", "no-such-model-avr");
      try {
        CommonFlags.getAvrModel(flags, false);
        fail();
      } catch (final InvalidParamsException e) {
        // expected
      }
    } finally {
      Diagnostic.setLogStream();
    }
  }

  public void testAvrModelDirFailures() throws Exception {
    final MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());
    try (TestDirectory tmpDir = new TestDirectory()) {
      try {
        final CFlags flags = new CFlags("blah", TestUtils.getNullPrintStream(), mps.printStream());
        final String oldAvrModels = System.getProperty(CommonFlags.ENVIRONMENT_MODELS_DIR);
        try {
          final File f = new File(tmpDir, "blah");
          try {
            System.setProperty(CommonFlags.ENVIRONMENT_MODELS_DIR, f.getPath());
            CommonFlags.getAvrModel(flags, false);
            fail();
          } catch (final InvalidParamsException e) {
            assertFalse(f.exists());
            assertEquals("The AVR models directory cannot be found or is not a directory: " + f.getPath(), e.getMessage());
          }
          assertTrue(f.createNewFile());
          try {
            System.setProperty(CommonFlags.ENVIRONMENT_MODELS_DIR, f.getPath());
            CommonFlags.getAvrModel(flags, false);
            fail();
          } catch (final InvalidParamsException e) {
            assertTrue(f.exists());
            assertEquals("The AVR models directory cannot be found or is not a directory: " + f.getPath(), e.getMessage());
          }
        } finally {
          if (oldAvrModels == null) {
            System.clearProperty(CommonFlags.ENVIRONMENT_MODELS_DIR);
          } else {
            System.setProperty(CommonFlags.ENVIRONMENT_MODELS_DIR, oldAvrModels);
          }
        }
      } finally {
        Diagnostic.setLogStream();
      }
    }
  }

}
