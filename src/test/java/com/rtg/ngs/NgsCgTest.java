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

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Collections;

import com.rtg.launcher.SequenceParams;
import com.rtg.mode.DnaUtils;
import com.rtg.mode.SequenceMode;
import com.rtg.reader.PrereadArm;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 * Tests of Ngs paired end capability.
 */
public class NgsCgTest extends TestCase {

  private File mDir;

  @Override
  public void setUp() throws Exception {
    Diagnostic.setLogStream();
    mDir = FileHelper.createTempDirectory();
  }

  @Override
  public void tearDown() {
    assertTrue(FileHelper.deleteAll(mDir));
    mDir = null;
  }

  /** the better getParams */
  NgsParams getParamsCG(OutputStream out, NgsMaskParams mask, NgsTestUtils.TestPairedEndParams test) throws IOException {
    final File subjectsDir = FileHelper.createTempDirectory(mDir);
    final File secondDir = FileHelper.createTempDirectory(mDir);
    final File queriesDir = FileHelper.createTempDirectory(mDir);
    final File hitsDir = FileHelper.createTempDirectory(mDir);
    //System.err.println("hitsDir=" + hitsDir);
    ReaderTestUtils.getReaderDNAFastqCG(test.mSubjects, subjectsDir, PrereadArm.LEFT).close();
    ReaderTestUtils.getReaderDNAFastqCG(test.mSubjectsSecond, secondDir, PrereadArm.RIGHT).close();
    ReaderTestUtils.getReaderDNA(test.mQueries, queriesDir, null).close();
    final SequenceParams subjectParams = SequenceParams.builder().directory(subjectsDir).useMemReader(true).mode(SequenceMode.UNIDIRECTIONAL).create();
    final SequenceParams secondParams = SequenceParams.builder().directory(secondDir).useMemReader(true).mode(SequenceMode.UNIDIRECTIONAL).create();
    final SequenceParams queryParams = SequenceParams.builder().directory(queriesDir).useMemReader(true).loadNames(true).create();
    final NgsFilterParams filterParams = NgsFilterParams.builder().outputFilter(OutputFilter.PAIRED_END).topN(0).errorLimit(10).create();

    final NgsOutputParams outputParams = new NgsTestUtils.OverriddenNgsOutputParams(NgsTestUtils.OverriddenNgsOutputParams.builder().outStream(out).progress(test.mProgress).outputDir(new File(hitsDir, "log")).filterParams(filterParams));
    return NgsParams.builder()
        .stepSize(test.mStepSize)
        .buildFirstParams(subjectParams)
        .buildSecondParams(secondParams)
        .searchParams(queryParams)
        .outputParams(outputParams)
        .maskParams(mask)
        .listeners(Collections.singleton(test.mListener))
        .maxFragmentLength(test.mMaxInsertSize)
        .minFragmentLength(test.mMinInsertSize)
        .create();
  }

  /** the better check */
  void checkCG(final NgsMaskParams mask, final NgsTestUtils.TestPairedEndParams test) throws Exception {
    final String actual;
    try (ByteArrayOutputStream out = new ByteArrayOutputStream()) {
      final NgsParams params = getParamsCG(out, mask, test);
      NgsTestUtils.execNgs(params);
      actual = TestUtils.stripSAMHeader(FileUtils.fileToString(new File(params.directory(), NgsOutputParams.MATED_SAM_FILE_NAME)));
    }
    //System.err.println("actual\n" + actual);
    assertTrue(TestUtils.sameLines(test.mExpected, actual, false));
  }

  private static final String SEQ0_CG0 = "TGACTGCA" + "TGCATGCA" + "CATTGCTA" + "CATTGCTA" + "TGC";

  private static final String SEQ1_CG0 = "GCCTAAGG" + "GCCTAAGG" + "GATACAAA" + "GATACAAA" + "GCC";

  private static final String QUAL_CG0 = "34563456" + "34563456" + "34563456" + "34563456" + "345";

  private static final String READ_FIRST_CG0 = ""
      + "@test" + LS
      + SEQ0_CG0 + LS
      + "+" + LS
      + QUAL_CG0 + LS
      ;

  private static final String READ_SECOND_CG0 = ""
      + "@test" + LS
      + SEQ1_CG0 + LS
      + "+" + LS
      + QUAL_CG0 + LS
      ;

  private static final String SPACED0_CG0 = "TGACTGCA" + "TGCATGCA" + "CATTGCTA" + "C" + "aaaaaa" + "ATTGCTA" + "TGC";

  private static final String SPACED1_CG0 = "GCCTAAGG" + "GC" + "aaaaaa" + "CTAAGG" + "GATACAAA" + "GATACAAA" + "GCC";

  private static final String GENOME_CG0 = ""
      + ">genome" + LS
      + SPACED0_CG0
      + "CCATTACGTGC" //spacer
      + SPACED1_CG0 + LS
      ;

  //One correctly mated pair - end to end test using CGMaska0.
  public void testCG0() throws Exception {
    Diagnostic.setLogStream();
    checkCG(new NgsMaskParamsExplicit("CGMaska0"),
        new NgsTestUtils.TestPairedEndParams(READ_FIRST_CG0, READ_SECOND_CG0, GENOME_CG0, FileHelper.resourceToString("com/rtg/ngs/resources/ngscg_0.txt"),
            0, 100, 0L)
        );
  }

  private static final String GENOME_CG1 = ""
      + ">genome" + LS
      + DnaUtils.reverseComplement(SPACED1_CG0)
      + "CCATTACGTGC" //spacer
      + DnaUtils.reverseComplement(SPACED0_CG0) + LS
      ;

  //One correctly mated pair - end to end test using CGMaska0.
  public void testCG1() throws Exception {
    Diagnostic.setLogStream();
    checkCG(new NgsMaskParamsExplicit("CGMaska0"),
        new NgsTestUtils.TestPairedEndParams(READ_FIRST_CG0, READ_SECOND_CG0, GENOME_CG1, FileHelper.resourceToString("com/rtg/ngs/resources/ngscg_1.txt"),
            0, 100, 0L)
        );
  }

  //carefully place two (aa) substitutions at the end of the 17-nt sub-window
  private static final String SEQ0_CG2 = "TGACTGCA" + "TGCATGCA" + "CA" + "AA" + "GCTA" + "CATTGCTA" + "TGC";

  private static final String SEQ1_CG2 = "GCCTAAGG" + "GCCTAAG" + "AA" + "ATACAAA" + "GATACAAA" + "GCC";

  private static final String READ_FIRST_CG2 = ""
      + "@test" + LS
      + SEQ0_CG2 + LS
      + "+" + LS
      + QUAL_CG0 + LS
      ;

  private static final String READ_SECOND_CG2 = ""
      + "@test" + LS
      + SEQ1_CG2 + LS
      + "+" + LS
      + QUAL_CG0 + LS
      ;

  //One correctly mated pair - end to end test using CGMaska0.
  public void testCG2() throws Exception {
    Diagnostic.setLogStream();
    checkCG(new NgsMaskParamsExplicit("CGMaska0"),
        new NgsTestUtils.TestPairedEndParams(READ_FIRST_CG2, READ_SECOND_CG2, GENOME_CG0, FileHelper.resourceToString("com/rtg/ngs/resources/ngscg_2.txt"),
            0, 100, 0L)
        );
  }
}
