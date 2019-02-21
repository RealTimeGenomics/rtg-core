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

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

import com.rtg.launcher.CommonFlags;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.sam.SamCommandHelper;
import com.rtg.util.IntegerOrPercentage;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.TestUtils;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.TestCFlags;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.test.FileHelper;
import com.rtg.util.test.MockEventListener;

import junit.framework.TestCase;

/**
 * Test class
 */
public class MapFlagsTest extends TestCase {

  private static final String LONG_READ_1 = "CATAATGACGGCTGGGCTACTGGACATATGTACGCGGTCTCCGGGAACAAGCAGCATGGAGTTTCCCCTGACGTATCGGTGATGTGATTGACATAACGATCAGATTTCAAAAGGAGTTCGCGCATTCCAGAGGACGCTATGCACGTTGGT";
  private static final String LONG_READ_2 = "ACAAGTACAAGGAGGCGGACTATCGCACTCCAAGTTAGCCGGTTTGAACATAGAAGCTCTCCCGAGCTCGCGCTATCCATCGATGTGGATGAGGAGTCAGCAATAAACTCGCTGATAGCCTGATCCAGTCAACTTTATCCATGTCAGGCA";

  private static final String TEMPLATE_LONG = LONG_READ_1 + "atcgaggtcatctagcagcatcatcgacttatcgacatctacgatcgagcgactatcgactatg" + LONG_READ_2;

  private CFlags mFlags = null;

  @Override
  public void setUp() {
    mFlags = new CFlags();
    MapFlags.initMaskFlagsOnly(mFlags);
    MapFlags.initPairedEndFlags(mFlags);
    MapFlags.initWordSize(mFlags, "my description");
    MapFlags.initMapFlags(mFlags);
    MapFlags.initReadFreqFlag(mFlags, 5);
    MapFlags.initSamOutputFlag(mFlags);
    SamCommandHelper.initSamRg(mFlags);
  }

  @Override
  public void tearDown() {
    mFlags = null;
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
      MapFlags.validateMaskParams(5, 6, 0, 0, 1);
      fail();
    } catch (final InvalidParamsException ipe) {
      assertEquals(ErrorType.WORD_NOT_LESS_READ, ipe.getErrorType());
//      ps.flush();
//      assertTrue(baos.toString().contains("The word length \"6\" should be less than the read length \"5\""));
    } finally {
      Diagnostic.setLogStream();
    }

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
      MapFlags.validateMaskParams(readlen, w, subs, indels, indellen);
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
    } finally {
      Diagnostic.setLogStream();
    }
  }

  public void testNegativeRepeatFreq() {
    final CFlags flags = new CFlags();
    CommonFlags.initRepeatFrequencyFlag(flags);
    flags.setFlags("-r", "-1");
    assertFalse(MapFlags.checkRepeatFrequency(flags));
    TestUtils.containsAll(flags.getParseMessage(), "--" + CommonFlags.REPEAT_FREQUENCY_FLAG + " must be in the range [1");

    flags.setFlags("-r", "0");
    assertFalse(MapFlags.checkRepeatFrequency(flags));
    TestUtils.containsAll(flags.getParseMessage(), "--" + CommonFlags.REPEAT_FREQUENCY_FLAG + " must be in the range [1");
  }


  public void testPercentRepeatFreq() {
    final MockEventListener ev = new MockEventListener();
    Diagnostic.addListener(ev);
    final CFlags flags = new CFlags();
    flags.registerOptional('r', CommonFlags.REPEAT_FREQUENCY_FLAG, IntegerOrPercentage.class, CommonFlags.INT, "maximum repeat frequency", new IntegerOrPercentage(20));
    try {
      flags.setFlags("-r", "-1");
      assertFalse(MapFlags.checkPercentRepeatFrequency(flags));
      assertTrue(ev.compareErrorMessage("Error: The specified flag \"--" + CommonFlags.REPEAT_FREQUENCY_FLAG + "\" has invalid value \"" + -1 + "\". It should be greater than or equal to \"1\"."));

      flags.setFlags("-r", "0");
      assertFalse(MapFlags.checkPercentRepeatFrequency(flags));
      assertTrue(ev.compareErrorMessage("Error: The specified flag \"--" + CommonFlags.REPEAT_FREQUENCY_FLAG + "\" has invalid value \"" + 0 + "\". It should be greater than or equal to \"1\"."));

      flags.setFlags("-r", "100001");
      assertFalse(MapFlags.checkPercentRepeatFrequency(flags));
      assertTrue(ev.compareErrorMessage("Error: The specified flag \"--" + CommonFlags.REPEAT_FREQUENCY_FLAG + "\" has invalid value \"" + 100001 + "\". It should be less than or equal to \"100000\"."));

      flags.setFlags("-r", "-1%");
      assertFalse(MapFlags.checkPercentRepeatFrequency(flags));
      assertTrue(ev.compareErrorMessage("Error: The specified flag \"--" + CommonFlags.REPEAT_FREQUENCY_FLAG + "\" has invalid value \"" + -1 + "%\". It should be greater than or equal to \"0%\"."));

      flags.setFlags("-r", "101%");
      assertFalse(MapFlags.checkPercentRepeatFrequency(flags));
      assertTrue(ev.compareErrorMessage("Error: The specified flag \"--" + CommonFlags.REPEAT_FREQUENCY_FLAG + "\" has invalid value \"" + 101 + "%\". It should be less than or equal to \"100%\"."));

      flags.setFlags("-r", "1");
      assertTrue(MapFlags.checkPercentRepeatFrequency(flags));
    } finally {
      Diagnostic.removeListener(ev);
    }
  }

  public void testMapIOFlags() {
    final CFlags flags = new CFlags("blah", TestUtils.getNullPrintStream(), TestUtils.getNullPrintStream());
    MapFlags.initMapIOFlags(flags);
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
    MapFlags.initMapIOFlags(flags);
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

      final File newDir = new File(tempDir, "blsdkjrlks");
      flags.setFlags("-o", newDir.getAbsolutePath());
      assertTrue(CommonFlags.validateOutputDirectory(flags));
      flags.setFlags("-o", tempDir.getAbsolutePath());
      assertFalse(CommonFlags.validateOutputDirectory(flags));
    } finally {
      Diagnostic.removeListener(ev);
      assertTrue(FileHelper.deleteAll(tempDir));
    }
  }

  public void testWordSize() {
    final CFlags flags = new CFlags();
    MapFlags.initWordSize(flags, "");
    assertEquals(5, MapFlags.getWordSize(flags, 10, 7));
    assertEquals(7, MapFlags.getWordSize(flags, 15, 7));
    flags.setFlags("--" + MapFlags.WORDSIZE_FLAG, "8");
    assertEquals(8, MapFlags.getWordSize(flags, 15, 7));
  }

  public void testThreads() {
    Diagnostic.setLogStream(TestUtils.getNullPrintStream());
    try {
      final CFlags flags = new CFlags("blah", TestUtils.getNullPrintStream(), TestUtils.getNullPrintStream());
      MapFlags.initSharedFlagsOnly(flags);
      flags.setFlags("-T", "0");
      assertFalse(CommonFlags.validateThreads(flags));
      flags.setFlags("-T", "-1");
      assertFalse(CommonFlags.validateThreads(flags));
      flags.setFlags("-T", Integer.toString(Integer.MAX_VALUE));
      assertFalse(CommonFlags.validateThreads(flags));
      flags.setFlags("-T", "10");
      assertTrue(CommonFlags.validateThreads(flags));
    } finally {
      Diagnostic.setLogStream();
    }
  }



  public void testInitMaskFlag() {
    final CFlags flags = new CFlags("blah", TestUtils.getNullPrintStream(), TestUtils.getNullPrintStream());
    MapFlags.initMaskFlags(flags);

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
    MapFlags.initSharedFlags(flags);

    assertNotNull(flags.getFlag(MapFlags.STEP_FLAG));
    assertNotNull(flags.getFlag(CommonFlags.THREADS_FLAG));
    assertNotNull(flags.getFlag(CommonFlags.NO_GZIP));

    TestUtils.containsAll(flags.getUsageString(),
      "-s", "step size (Default is word size)",
      "maximum repeat frequency",
      "-T", "number of threads",
      "-Z", "do not gzip the output"
    );
  }

  public void testCheckFlags() {
    TestCFlags.check(mFlags);
    TestCFlags.checkUsage(mFlags,
       "guaranteed minimum number of substitutions which will be detected (Default is 1)",
       "guaranteed minimum number of indels which will be detected (Default is 1)",
       "maximum mismatches for mappings in single-end mode (as absolute value or percentage of read length) (Default is 10%)",
       "maximum mismatches for mappings of unmated results (as absolute value or percentage of read length) (Default is 10%)",
       "maximum mismatches for mappings across mated results, alias for ",
       "--max-mismatches (as absolute value or percentage of read length) (Default is 10%)",
       "directory used for temporary files (Defaults to output directory)",
       "--legacy-cigars ", "use legacy cigars in output",
       "--read-names ", "use read name in output instead of read id (Uses more RAM)",
       "--sam ", "output the alignment files in SAM format",
       "--sam-rg", "file containing a single valid read group SAM header line");

    assertEquals(CommonFlagCategories.REPORTING, mFlags.getFlag(MapFlags.MAX_ALIGNMENT_MISMATCHES).getCategory());
    assertEquals(CommonFlagCategories.UTILITY, mFlags.getFlag(MapFlags.LEGACY_CIGARS).getCategory());
    assertEquals(CommonFlagCategories.SENSITIVITY_TUNING, mFlags.getFlag(MapFlags.READ_FREQUENCY_FLAG).getCategory());
  }

  public void testWord() {
    checkInRange("--word", 0, 1);
  }

  public void testSubs() throws Exception {
    all4combo("-a", -1, 0);
  }

  public void testIndel() throws Exception {
    all4combo("-b", -1, 0);
  }

  public void testIndelLength() throws Exception {
    all4combo("-c", -1, 1);
  }

  public void testMaxFragmentSize() {
    checkInRange("--max-fragment-size", 0, 1);
  }

  public void testMinFragmentSize() {
    checkInRange("--min-fragment-size", -1, 0);
  }

  public void testMaxMatedAlignmentScore() throws Exception {
    all4combo("--max-mated-mismatches", -1, 0);
  }

  public void testMaxAlignmentScore() throws Exception {
    all4combo("--max-mismatches", -1, 0);
  }

  public void testMaxUnmatedAlignmentScore() throws Exception {
    all4combo("--max-unmated-mismatches", -1, 0);
  }

  private void all4combo(final String flag, final int value, final int minValue) {
    checkMapParams(new String[] {flag, "" + value}, new String[] {"The specified flag \"" + flag + "\" has invalid value \"" + value + "\". It should be greater than or equal to \"" + minValue + "\"."});
  }

  private void checkInRange(final String flag, final int value, final int minValue) {
    mFlags.setFlags(flag, "" + value);
    MapFlags.validateMapParams(mFlags);
    TestUtils.containsAll(mFlags.getParseMessage(), flag + " must be in the range [" + minValue);
  }


  private void checkMapParams(final String[] additionalArgs, final String[] diagMessage) {
    final MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());
    try {
      mFlags.setFlags(additionalArgs);
      MapFlags.validateMapParams(mFlags);
      TestUtils.containsAll(mps.toString(), diagMessage);
    } finally {
      Diagnostic.setLogStream();
    }
  }
}
