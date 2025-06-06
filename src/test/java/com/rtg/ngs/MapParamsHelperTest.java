/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
package com.rtg.ngs;

import static com.rtg.reader.SamBamSequenceDataSourceTest.SAM_NL;

import java.io.File;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.SequenceParams;
import com.rtg.mode.SequenceMode;
import com.rtg.ngs.MapParamsHelper.FastaSequenceParamsCallable;
import com.rtg.ngs.MapParamsHelper.SamSequenceParamsCallable;
import com.rtg.ngs.MapParamsHelper.SdfSequenceParamsCallable;
import com.rtg.ngs.NgsFilterParams.NgsFilterParamsBuilder;
import com.rtg.reader.DataSourceDescription;
import com.rtg.reader.PrereadArm;
import com.rtg.reader.QualityFormat;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SamBamSequenceDataSourceTest;
import com.rtg.reader.SdfId;
import com.rtg.reader.SimpleNames;
import com.rtg.reader.SourceFormat;
import com.rtg.reference.Sex;
import com.rtg.sam.SamCommandHelper;
import com.rtg.util.IntegerOrPercentage;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;

import htsjdk.samtools.SAMReadGroupRecord;

/**
 */
public class MapParamsHelperTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new MapCli();
  }

  public void testCommonParams() throws Exception {
    final MemoryPrintStream err = new MemoryPrintStream();
    Diagnostic.setLogStream(err.printStream());
    final NgsParamsBuilder npb = new NgsParamsBuilder().stepSize(3);
    final CFlags flags = new CFlags();
    try (final TestDirectory tmpDir = new TestDirectory()) {
      MapFlags.initInputFormatFlags(flags);
      MapFlags.initMapIOFlags(flags);
      MapFlags.initSharedFlagsOnly(flags);
      CommonFlags.initReadRange(flags);
      MapFlags.initStepSize(flags);
      flags.registerOptional(MapFlags.THREAD_MULTIPLIER, Integer.class, CommonFlags.INT, "").setCategory(CommonFlagCategories.UTILITY);
      MapFlags.initReadFreqFlag(flags, 3);

      ReaderTestUtils.createPairedReaderDNA(">s" + StringUtils.LS + "AGT", ">s" + StringUtils.LS + "TGA", tmpDir, new SdfId(5L));

      final File r = new File(tmpDir, "right");
      flags.setFlags("-i", tmpDir.getPath(),
                            "-t", r.getPath(),
                            "-o", tmpDir.getPath(),
                            "--" + MapFlags.THREAD_MULTIPLIER, "3",
                            "--" + CommonFlags.REPEAT_FREQUENCY_FLAG, "33%",
                            "--" + CommonFlags.END_READ_ID, "2");

      err.reset();
      assertEquals(3, MapParamsHelper.populateCommonMapParams(npb, flags, 12, 3, false, false));
      assertEquals(2, err.toString().split("The end sequence id \"2\" is out of range, it must be from \"1\" to \"1\". Defaulting end to \"1\"").length - 1);

      flags.setFlags("-i", tmpDir.getPath(),
                            "-t", r.getPath(),
                            "-o", tmpDir.getPath(),
                            "--" + MapFlags.THREAD_MULTIPLIER, "1",
                            "--" + CommonFlags.THREADS_FLAG, "1",
                            "--" + CommonFlags.REPEAT_FREQUENCY_FLAG, "33%",
                            "--" + CommonFlags.END_READ_ID, "3");
      assertEquals(3, MapParamsHelper.populateCommonMapParams(npb, flags, 8, 2, false, false));
      assertEquals(2, err.toString().split("The end sequence id \"3\" is out of range, it must be from \"1\" to \"1\". Defaulting end to \"1\"").length - 1);
      assertTrue(npb.mBuildFirstParams.readerParams().toString().contains("usememory=" + true));
      assertTrue(npb.mBuildSecondParams.readerParams().toString().contains("usememory=" + true));

      try {
        npb.stepSize(4);
        MapParamsHelper.populateCommonMapParams(npb, flags, 8, 2, false, false);
        fail();
      } catch (InvalidParamsException ipe) {
        assertEquals("Step size (4) must be less than or equal to max read length (3)", ipe.getMessage());
      }
    }
  }

  public void testLongCGErrorMessage() throws Exception {
    try (final TestDirectory mainOut = new TestDirectory()) {
      final File left = new File(mainOut, "left");
      final String inputDnaSequence = "@test\nacgtacgtacgtacgtacgtacgtacgtacgtacg\n+\n###################################";
      ReaderTestUtils.getReaderDNAFastqCG(inputDnaSequence, left, PrereadArm.LEFT);
      final File right = new File(mainOut, "right");
      ReaderTestUtils.getReaderDNAFastqCG(inputDnaSequence, right, PrereadArm.RIGHT);
      final File template = new File(mainOut, "template");
      ReaderTestUtils.getReaderDNA(">template\nacgt", template, null);
      final File out = new File(mainOut, "out");
      final MapCli map = (MapCli) mCli;

      checkHandleFlagsOut("-i", mainOut.getPath(),
          "-t", template.getPath(),
          "-o", out.getPath());

      final MemoryPrintStream mps = new MemoryPrintStream();
      Diagnostic.setLogStream(mps.printStream());

      try {
        map.makeParams();
        fail();
      } catch (final InvalidParamsException ipe) {
        assertEquals(ErrorType.IS_A_CG_SDF, ipe.getErrorType());
//    assertTrue(mps.toString().contains("The file \"" + mainOut.getPath() + "\" contains Complete Genomics reads, please use \"cgmap\" module to map these reads."));
      }
    }
  }

  public void testMaskParams() {
    final MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());
    try {
      final CFlags flags = new CFlags();
      MapFlags.initMaskFlagsOnly(flags);
      MapFlags.initWordSize(flags, "");
      try {
        flags.setFlags("-w", "37", "-a", "1", "-b", "1", "-c", "1");

        MapParamsHelper.makeMaskParams(flags, 35, false, 18);
        fail();
      } catch (final InvalidParamsException ipe) {
        assertEquals(ErrorType.WORD_NOT_LESS_READ, ipe.getErrorType());
//        assertTrue(mps.toString(), mps.toString().contains("The word length \"37\" should be less than the read length \"35\"."));
      }
      try {
        mps.reset();
        flags.setFlags("-w", "32", "-a", "1", "-b", "1", "-c", "5");

        MapParamsHelper.makeMaskParams(flags, 35, false, 32);
        fail();
      } catch (final InvalidParamsException ipe) {
        assertEquals(ErrorType.INVALID_MAX_INTEGER_FLAG_VALUE, ipe.getErrorType());
//        assertTrue(mps.toString(), mps.toString().contains("The specified flag \"-c\" has invalid value \"5\". It should be less than or equal to \"3\""));
      }
      try {
        mps.reset();
        flags.setFlags("-w", "32", "-a", "1", "-b", "1", "-c", "0");

        MapParamsHelper.makeMaskParams(flags, 35, false, 32);
        fail();
      } catch (final InvalidParamsException ipe) {
        assertEquals(ErrorType.INVALID_MIN_INTEGER_FLAG_VALUE, ipe.getErrorType());
//        assertTrue(mps.toString(), mps.toString().contains("The specified flag \"-c\" has invalid value \"0\". It should be greater than or equal to \"1\""));
      }
      try {
        mps.reset();
        flags.setFlags("-w", "32", "-a", "1", "-b", "5", "-c", "1");

        MapParamsHelper.makeMaskParams(flags, 35, false, 32);
        fail();
      } catch (final InvalidParamsException ipe) {
        assertEquals(ErrorType.INVALID_MAX_INTEGER_FLAG_VALUE, ipe.getErrorType());
//        assertTrue(mps.toString(), mps.toString().contains("The specified flag \"-b\" has invalid value \"5\". It should be less than or equal to \"3\""));
      }
      try {
        mps.reset();
        flags.setFlags("-w", "32", "-a", "5", "-b", "1", "-c", "1");

        MapParamsHelper.makeMaskParams(flags, 35, false, 32);
        fail();
      } catch (final InvalidParamsException ipe) {
        assertEquals(ErrorType.INVALID_MAX_INTEGER_FLAG_VALUE, ipe.getErrorType());
//        assertTrue(mps.toString(), mps.toString().contains("The specified flag \"-a\" has invalid value \"5\". It should be less than or equal to \"3\""));
      }
      try {
        mps.reset();
        flags.setFlags("-w", "32", "-a", "1", "-b", "1", "-c", "1");

        MapParamsHelper.makeMaskParams(flags, 35, false, 18);
        fail();
      } catch (final InvalidParamsException ipe) {
        assertEquals(ErrorType.INVALID_MASK_PARAMS, ipe.getErrorType());
//        assertTrue(mps.toString(), mps.toString().contains("The combination of w, a, b, and c is not valid for this data set. Please consult the user manual."));
      }
      flags.setFlags("-w", "8", "-a", "1", "-b", "1", "-c", "1");
      assertNotNull(MapParamsHelper.makeMaskParams(flags, 35, false, 18));
    } finally {
      Diagnostic.setLogStream();
    }
  }



  public void testPopulateAlignmentScores() {
    final CFlags flags = new CFlags();
    MapFlags.initMaskFlagsOnly(flags);
    MapFlags.initMapFlags(flags);
    final NgsFilterParamsBuilder npb = new NgsFilterParamsBuilder();

    flags.setFlags("--" + MapFlags.MAX_ALIGNMENT_MISMATCHES, "3");

    MapParamsHelper.populateAlignmentScoreSettings(flags, npb, false, null);

    assertEquals(new IntegerOrPercentage(3), npb.mMatedMaxMismatches);
    assertEquals(new IntegerOrPercentage(3), npb.mUnmatedMaxMismatches);

    flags.setFlags("--" + MapFlags.MAX_ALIGNMENT_MISMATCHES, "5%",
                          "--" + MapFlags.UNMATED_MISMATCH_THRESHOLD, "7%");
    final MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());
    try {
      MapParamsHelper.populateAlignmentScoreSettings(flags, npb, true, null);

      assertEquals(new IntegerOrPercentage("5%"), npb.mMatedMaxMismatches);
      assertEquals(new IntegerOrPercentage("7%"), npb.mUnmatedMaxMismatches);

      TestUtils.containsAll(mps.toString(), "--" + MapFlags.UNMATED_MISMATCH_THRESHOLD + " should not be greater than --" + MapFlags.MATED_MISMATCH_THRESHOLD);
    } finally {
      Diagnostic.setLogStream();
    }
  }

  public void testPopulateAlignmentScoresDefaults() {
    final CFlags flags = new CFlags();

    flags.registerOptional('e', MapFlags.MAX_ALIGNMENT_MISMATCHES, IntegerOrPercentage.class, CommonFlags.INT, "maximum alignment score for mappings in single-end mode (as absolute value or percentage of read length)", IntegerOrPercentage.valueOf("5%")).setCategory(CommonFlagCategories.REPORTING);
    flags.registerOptional('E', MapFlags.UNMATED_MISMATCH_THRESHOLD, IntegerOrPercentage.class, CommonFlags.INT, "maximum alignment score for mappings of unmated results (as absolute value or percentage of read length)", IntegerOrPercentage.valueOf("10%")).setCategory(CommonFlagCategories.REPORTING);
    flags.registerOptional(MapFlags.MATED_MISMATCH_THRESHOLD, IntegerOrPercentage.class, CommonFlags.INT, "maximum alignment score for mappings across mated results, alias for --" + MapFlags.MAX_ALIGNMENT_MISMATCHES + " (as absolute value or percentage of read length)", IntegerOrPercentage.valueOf("15%")).setCategory(CommonFlagCategories.REPORTING);


    final NgsFilterParamsBuilder npb = new NgsFilterParamsBuilder();
    MapParamsHelper.populateAlignmentScoreSettings(flags, npb, false, null);

    assertEquals(new IntegerOrPercentage("5%"), npb.mMatedMaxMismatches);
    assertEquals(new IntegerOrPercentage("5%"), npb.mUnmatedMaxMismatches);

    MapParamsHelper.populateAlignmentScoreSettings(flags, npb, true, null);

    assertEquals(new IntegerOrPercentage("5%"), npb.mMatedMaxMismatches); //I think this behaviour is dumb, but so be it.
    assertEquals(new IntegerOrPercentage("10%"), npb.mUnmatedMaxMismatches);
  }

  public void testPopulateAlignmentScoresDefaultsIonTorrent() throws Exception {
    try (final TestDirectory outer = new TestDirectory()) {
        final String rgstring = "@RG\tID:L23\tSM:NA123\tPL:IONTORRENT\n";
        final File header = new File(outer, "header");
        FileUtils.stringToFile(rgstring, header);

        final CFlags flags = new CFlags();

        SamCommandHelper.initSamRg(flags);
        flags.registerOptional('e', MapFlags.MAX_ALIGNMENT_MISMATCHES, IntegerOrPercentage.class, CommonFlags.INT, "maximum alignment score for mappings in single-end mode (as absolute value or percentage of read length)", IntegerOrPercentage.valueOf("5%")).setCategory(CommonFlagCategories.REPORTING);
        flags.registerOptional('E', MapFlags.UNMATED_MISMATCH_THRESHOLD, IntegerOrPercentage.class, CommonFlags.INT, "maximum alignment score for mappings of unmated results (as absolute value or percentage of read length)", IntegerOrPercentage.valueOf("10%")).setCategory(CommonFlagCategories.REPORTING);
        flags.registerOptional(MapFlags.MATED_MISMATCH_THRESHOLD, IntegerOrPercentage.class, CommonFlags.INT, "maximum alignment score for mappings across mated results, alias for --" + MapFlags.MAX_ALIGNMENT_MISMATCHES + " (as absolute value or percentage of read length)", IntegerOrPercentage.valueOf("15%")).setCategory(CommonFlagCategories.REPORTING);

        flags.setFlags("--" + SamCommandHelper.SAM_RG, header.getPath());

        final NgsFilterParamsBuilder npb = new NgsFilterParamsBuilder();
        final SAMReadGroupRecord rg = SamCommandHelper.validateAndCreateSamRG(rgstring, SamCommandHelper.ReadGroupStrictness.REQUIRED);
        MapParamsHelper.populateAlignmentScoreSettings(flags, npb, false, rg);

        assertEquals(new IntegerOrPercentage("10%"), npb.mMatedMaxMismatches);
        assertEquals(new IntegerOrPercentage("10%"), npb.mUnmatedMaxMismatches);

        MapParamsHelper.populateAlignmentScoreSettings(flags, npb, true, rg);

        assertEquals(new IntegerOrPercentage("10%"), npb.mMatedMaxMismatches); //I think this behaviour is dumb, but so be it.
        assertEquals(new IntegerOrPercentage("10%"), npb.mUnmatedMaxMismatches);
    }
  }
  public void testPopulateAlignmentScoresIonTorrent() throws Exception {
    try (final TestDirectory outer = new TestDirectory()) {
      final String rgstring = "@RG\tID:L23\tSM:NA123\tPL:IONTORRENT\n";
      final File header = new File(outer, "header");
      FileUtils.stringToFile(rgstring, header);

      final CFlags flags = new CFlags();

      SamCommandHelper.initSamRg(flags);
      flags.registerOptional('e', MapFlags.MAX_ALIGNMENT_MISMATCHES, IntegerOrPercentage.class, CommonFlags.INT, "maximum alignment score for mappings in single-end mode (as absolute value or percentage of read length)", IntegerOrPercentage.valueOf("5%")).setCategory(CommonFlagCategories.REPORTING);
      flags.registerOptional('E', MapFlags.UNMATED_MISMATCH_THRESHOLD, IntegerOrPercentage.class, CommonFlags.INT, "maximum alignment score for mappings of unmated results (as absolute value or percentage of read length)", IntegerOrPercentage.valueOf("10%")).setCategory(CommonFlagCategories.REPORTING);
      flags.registerOptional(MapFlags.MATED_MISMATCH_THRESHOLD, IntegerOrPercentage.class, CommonFlags.INT, "maximum alignment score for mappings across mated results, alias for --" + MapFlags.MAX_ALIGNMENT_MISMATCHES + " (as absolute value or percentage of read length)", IntegerOrPercentage.valueOf("15%")).setCategory(CommonFlagCategories.REPORTING);

      flags.setFlags("--" + SamCommandHelper.SAM_RG, header.getPath(), "-e", "2%", "-E", "3%", "--" + MapFlags.MATED_MISMATCH_THRESHOLD, "4%");

      final NgsFilterParamsBuilder npb = new NgsFilterParamsBuilder();
      final SAMReadGroupRecord rg = SamCommandHelper.validateAndCreateSamRG(rgstring, SamCommandHelper.ReadGroupStrictness.REQUIRED);
      MapParamsHelper.populateAlignmentScoreSettings(flags, npb, false, rg);

      assertEquals(new IntegerOrPercentage("2%"), npb.mMatedMaxMismatches);
      assertEquals(new IntegerOrPercentage("2%"), npb.mUnmatedMaxMismatches);

      MapParamsHelper.populateAlignmentScoreSettings(flags, npb, true, rg);

      assertEquals(new IntegerOrPercentage("4%"), npb.mMatedMaxMismatches); //I think this behaviour is dumb, but so be it.
      assertEquals(new IntegerOrPercentage("3%"), npb.mUnmatedMaxMismatches);
    }
  }

  public void testSequenceParamsCallableSdf() throws Exception {
    try (final TestDirectory tmpDir = new TestDirectory()) {
      final MemoryPrintStream mps = new MemoryPrintStream();
      Diagnostic.setLogStream(mps.printStream());
      final File template = new File(tmpDir, "template");
      ReaderTestUtils.getReaderDNA(">template\nacgt", template, null);

      final SdfSequenceParamsCallable spc = new SdfSequenceParamsCallable(template, LongRange.NONE, new NameParams(false, false), SequenceMode.UNIDIRECTIONAL);
      SequenceParams sp = null;
      try {
        sp = spc.call();
        assertEquals(SequenceMode.UNIDIRECTIONAL, sp.mode());
        assertTrue(sp.readerParams().toString().contains("usememory=" + true));
        assertFalse(mps.toString().contains("Sequence names passed checksum"));
      } finally {
        if (sp != null) {
          sp.close();
        }
      }

      mps.reset();
      sp = null;
      final SdfSequenceParamsCallable spc2 = new SdfSequenceParamsCallable(template, LongRange.NONE, false, Sex.MALE, false, SequenceMode.BIDIRECTIONAL);
      try {
        sp = spc2.call();
        assertEquals(SequenceMode.BIDIRECTIONAL, sp.mode());
        assertTrue(sp.readerParams().toString().contains("usememory=" + false));
        assertTrue(mps.toString().contains("Sequence names passed checksum"));
        assertEquals(Sex.MALE, sp.sex());
      } finally {
        if (sp != null) {
          sp.close();
        }
      }

      mps.reset();
      sp = null;
      final SdfSequenceParamsCallable spc3 = new SdfSequenceParamsCallable(template, LongRange.NONE,  new NameParams(true, false), SequenceMode.BIDIRECTIONAL);
      try {
        sp = spc3.call();
        assertEquals(SequenceMode.BIDIRECTIONAL, sp.mode());
        assertTrue(sp.readerParams().toString().contains("usememory=" + true));
        assertTrue(mps.toString().contains("Sequence names passed checksum"));
        assertNull(sp.sex());
      } finally {
        if (sp != null) {
          sp.close();
        }
      }
    }
  }

  private static final DataSourceDescription FASTQ_DS = new DataSourceDescription(SourceFormat.FASTQ, QualityFormat.SANGER, false, false, false);

  public void testSequenceParamsCallableFasta() throws Exception {
    try (final TestDirectory tmpDir = new TestDirectory()) {
      final MemoryPrintStream mps = new MemoryPrintStream();
      Diagnostic.setLogStream(mps.printStream());
      final String inputDnaSequence = "@test\nacgtacgtacgtacgtacgtacgtacgtacgtacg\n+\n###################################\n";
      final File single = new File(tmpDir, "single.fastq");
      FileUtils.stringToFile(inputDnaSequence, single);

      final FastaSequenceParamsCallable spc = new FastaSequenceParamsCallable(single, FASTQ_DS, LongRange.NONE, PrereadArm.LEFT,  null, null, true, SequenceMode.UNIDIRECTIONAL);
      final SequenceParams sp = spc.call()[0];
      assertNull(spc.call()[1]);
      assertEquals(SequenceMode.UNIDIRECTIONAL, sp.mode());
      assertTrue(sp.readerParams().toString().contains("usememory=" + true));
      assertFalse(mps.toString().contains("Sequence names passed checksum"));
      TestUtils.containsAll(mps.toString(), "Processing left arm \"" + single.getPath());
      assertEquals(PrereadArm.UNKNOWN, sp.reader().getArm());
      assertTrue(sp.reader().getSdfId().equals(new SdfId(0L)));
      assertTrue(sp.reader().numberSequences() > 0);
      assertTrue(sp.reader().hasQualityData());
      assertEquals(single, sp.reader().path());

      mps.reset();
      final FastaSequenceParamsCallable spc2 = new FastaSequenceParamsCallable(single, FASTQ_DS, LongRange.NONE, PrereadArm.RIGHT,  null, null, true, SequenceMode.UNIDIRECTIONAL);
      final SequenceParams sp2 = spc2.call()[0];
      assertNull(spc2.call()[1]);
      assertEquals(SequenceMode.UNIDIRECTIONAL, sp2.mode());
      assertFalse(mps.toString().contains("Sequence names passed checksum"));
      assertTrue(mps.toString().contains("Processing right arm \"" + single.getPath()));
      assertEquals(PrereadArm.UNKNOWN, sp2.reader().getArm());
      assertEquals(single, sp2.reader().path());

      final File both = new File(tmpDir, "alternating.fastq");
      FileUtils.stringToFile(inputDnaSequence + inputDnaSequence, both);
      mps.reset();
      final DataSourceDescription fastqIntDs = new DataSourceDescription(SourceFormat.FASTQ, QualityFormat.SANGER, true, true, false);
      final FastaSequenceParamsCallable spc3 = new FastaSequenceParamsCallable(both, fastqIntDs, LongRange.NONE, PrereadArm.UNKNOWN,  null, null, true, SequenceMode.UNIDIRECTIONAL);
      final SequenceParams sp3l = spc3.call()[0];
      final SequenceParams sp3r = spc3.call()[1];
      assertEquals(SequenceMode.UNIDIRECTIONAL, sp3l.mode());
      assertEquals(SequenceMode.UNIDIRECTIONAL, sp3r.mode());
      assertFalse(mps.toString().contains("Sequence names passed checksum"));
      assertEquals(PrereadArm.LEFT, sp3l.reader().getArm());
      assertEquals(PrereadArm.RIGHT, sp3r.reader().getArm());
      assertEquals(both, sp3l.reader().path());
      assertEquals(both, sp3r.reader().path());
    }
  }

  public void testSequenceParamsCallableSam() throws Exception {
    try (final TestDirectory tmpDir = new TestDirectory()) {
      final MemoryPrintStream mps = new MemoryPrintStream();
      Diagnostic.setLogStream(mps.printStream());
      final File input = new File(tmpDir, "input");
      String inputSequence = SamBamSequenceDataSourceTest.SAM_HEADER + String.format(SamBamSequenceDataSourceTest.SAM_LINE_SINGLE, SAM_NL, "read0", "ACTG", "5555");
      FileUtils.stringToFile(inputSequence, input);

      final SamSequenceParamsCallable sps = new SamSequenceParamsCallable(input, new DataSourceDescription(SourceFormat.SAM, QualityFormat.SANGER, false, false, false), LongRange.NONE,  null, new SimpleNames(), true, SequenceMode.UNIDIRECTIONAL, new SamSequenceReaderParams(false, false));
      final SequenceParams[] sp = sps.call();
      assertEquals(2, sp.length);
      assertNull(sp[1]);
      assertNotNull(sp[0]);
      assertEquals(SequenceMode.UNIDIRECTIONAL, sp[0].mode());
      assertTrue(sp[0].readerParams().toString().contains("usememory=" + true));
      assertFalse(mps.toString().contains("Sequence names passed checksum"));
      assertEquals(PrereadArm.UNKNOWN, sp[0].reader().getArm());
      assertTrue(sp[0].reader().getSdfId().equals(new SdfId(0L)));
      assertTrue(sp[0].reader().numberSequences() > 0);
      assertTrue(sp[0].reader().hasQualityData());
      assertEquals(input, sp[0].reader().path());

      mps.reset();
      inputSequence = SamBamSequenceDataSourceTest.SAM_HEADER + String.format(SamBamSequenceDataSourceTest.SAM_LINE_LEFT, SAM_NL, "read0", "ACTG", "5555")
                                                              + String.format(SamBamSequenceDataSourceTest.SAM_LINE_RIGHT, SAM_NL, "read0", "GTCA", "6666");
      FileUtils.stringToFile(inputSequence, input);
      final SamSequenceParamsCallable sps2 = new SamSequenceParamsCallable(input, new DataSourceDescription(SourceFormat.SAM, QualityFormat.SANGER, false, true, false), LongRange.NONE,  null, new SimpleNames(), true, SequenceMode.UNIDIRECTIONAL, new SamSequenceReaderParams(false, false));
      final SequenceParams[] sp2 = sps2.call();
      assertNotNull(sp2[0]);
      assertNotNull(sp2[1]);
      assertEquals(SequenceMode.UNIDIRECTIONAL, sp2[0].mode());
      assertEquals(SequenceMode.UNIDIRECTIONAL, sp2[1].mode());
      assertFalse(mps.toString().contains("Sequence names passed checksum"));
      assertEquals(PrereadArm.LEFT, sp2[0].reader().getArm());
      assertEquals(PrereadArm.RIGHT, sp2[1].reader().getArm());
      assertEquals(input, sp2[0].reader().path());
      assertEquals(input, sp2[1].reader().path());

    }
  }

  public void testAlignmentFlags() {
    final CFlags flags = new CFlags();
    MapFlags.initAlignerPenaltyFlags(flags);

    flags.setFlags("--" + MapFlags.SOFT_CLIP_DISTANCE_FLAG, "9",
                          "--" + MapFlags.MISMATCH_PENALTY_FLAG, "7",
                          "--" + MapFlags.GAP_OPEN_PENALTY_FLAG, "2",
                          "--" + MapFlags.GAP_EXTEND_PENALTY_FLAG, "3");


    final NgsParamsBuilder b = new NgsParamsBuilder();

    assertEquals(0.5, b.mAlignerBandWidthFactor.getFactor(), 0.001);
    assertEquals(0, b.mIndelSoftClipDistance);
    assertEquals(9, b.mSubstitutionPenalty);
    assertEquals(19, b.mGapOpenPenalty);
    assertEquals(1, b.mGapExtendPenalty);

    assertEquals(b, MapParamsHelper.populateAlignerPenaltiesParams(b, flags));

    assertEquals(0.5, b.mAlignerBandWidthFactor.getFactor(), 0.001);
    assertEquals(9, b.mIndelSoftClipDistance);
    assertEquals(7, b.mSubstitutionPenalty);
    assertEquals(2, b.mGapOpenPenalty);
    assertEquals(3, b.mGapExtendPenalty);

    flags.setFlags("--" + MapFlags.ALIGNER_BAND_WIDTH_FACTOR_FLAG, "0.66",
                          "--" + MapFlags.SOFT_CLIP_DISTANCE_FLAG, "9",
                          "--" + MapFlags.MISMATCH_PENALTY_FLAG, "7",
                          "--" + MapFlags.GAP_OPEN_PENALTY_FLAG, "2",
                          "--" + MapFlags.GAP_EXTEND_PENALTY_FLAG, "3");

    assertEquals(b, MapParamsHelper.populateAlignerPenaltiesParams(b, flags));

    assertNotNull(b.mAlignerBandWidthFactor);
    assertEquals(0.66, b.mAlignerBandWidthFactor.getFactor(), 0.001);
    assertEquals(9, b.mIndelSoftClipDistance);
    assertEquals(7, b.mSubstitutionPenalty);
    assertEquals(2, b.mGapOpenPenalty);
    assertEquals(3, b.mGapExtendPenalty);
  }
}
