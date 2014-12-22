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

import com.rtg.launcher.CommonFlags;
import com.rtg.sam.SamCommandHelper;
import com.rtg.util.TestUtils;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.TestCFlags;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.MemoryPrintStream;

import junit.framework.TestCase;

/**
 * Test class
 */
public class MapFlagsTest extends TestCase {

  private CFlags mFlags = null;

  @Override
  public void setUp() {
    mFlags = new CFlags();
    CommonFlags.initMaskFlagsOnly(mFlags);
    CommonFlags.initPairedEndFlags(mFlags);
    CommonFlags.initWordSize(mFlags, "my description");
    MapFlags.initMapFlags(mFlags);
    MapFlags.initReadFreqFlag(mFlags, 5);
    MapFlags.initSamOutputFlag(mFlags);
    SamCommandHelper.initSamRg(mFlags);
  }

  @Override
  public void tearDown() {
    mFlags = null;
  }

  public void testCheckFlags() {
    TestCFlags.check(mFlags,
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

    assertEquals(CommonFlagCategories.REPORTING, mFlags.getFlag(CommonFlags.MAX_ALIGNMENT_MISMATCHES).getCategory());
    assertEquals(CommonFlagCategories.UTILITY, mFlags.getFlag(MapFlags.LEGACY_CIGARS).getCategory());
    assertEquals(CommonFlagCategories.SENSITIVITY_TUNING, mFlags.getFlag(MapFlags.READ_FREQUENCY_FLAG).getCategory());
  }

  public void testWord() throws Exception {
    all4combo("--word", 0, 1);
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

  public void testMaxFragmentSize() throws Exception {
    pairedEndCombo("--max-fragment-size", 0, 1);
  }

  public void testMinFragmentSize() throws Exception {
    pairedEndCombo("--min-fragment-size", -1, 0);
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

  private void all4combo(final String flag, final int value, final int minValue) throws Exception  {
    pairedEndCombo(flag, value, minValue);
  }

  private void pairedEndCombo(final String flag, final int value, final int minValue) {
    checkMapParams(new String[] {flag, "" + value}, new String[] {"The specified flag \"" + flag + "\" has invalid value \"" + value + "\". It should be greater than or equal to \"" + minValue + "\"."});
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
