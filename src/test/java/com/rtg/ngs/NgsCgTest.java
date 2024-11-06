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
