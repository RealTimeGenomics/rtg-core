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
import static com.rtg.util.StringUtils.TAB;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.Collections;

import com.rtg.launcher.SequenceParams;
import com.rtg.mode.SequenceMode;
import com.rtg.ngs.NgsTestUtils.OverriddenNgsOutputParams;
import com.rtg.ngs.NgsTestUtils.ParamsParams;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ListenerType;
import com.rtg.util.test.FileHelper;

/**
 */
public class NgsThreadTest extends NgsTaskFunctionalTest {

  @Override
  NgsParams getParams(final OutputStream out, final NgsMaskParams mask, final ParamsParams restParams, final ListenerType listener, final OutputFilter filter, final int topN, int numThreads) throws IOException {
    final File subjectsDir = FileHelper.createTempDirectory(mDir);
    final File queriesDir = FileHelper.createTempDirectory(mDir);
    final File hitsDir = FileHelper.createTempDirectory(mDir);
    //System.err.println("hitsDir=" + hitsDir);
    ReaderTestUtils.getReaderDNA(restParams.mSubjects, subjectsDir, null).close();
    ReaderTestUtils.getReaderDNA(restParams.mQueries, queriesDir, null).close();
    final SequenceParams subjectParams = SequenceParams.builder().directory(subjectsDir).mode(SequenceMode.UNIDIRECTIONAL).create();
    final SequenceParams queryParams = SequenceParams.builder().directory(queriesDir).create();
    final NgsFilterParams filterParams = NgsFilterParams.builder().outputFilter(filter).zip(restParams.mZip).topN(topN).errorLimit(restParams.mErrorLimit).create();
    final NgsOutputParams outputParams = new OverriddenNgsOutputParams(OverriddenNgsOutputParams.builder().outStream(out).progress(restParams.mProgress).outputDir(new File(hitsDir, "log")).filterParams(filterParams));
    return NgsParams.builder().numberThreads(4).buildFirstParams(subjectParams).searchParams(queryParams).outputParams(outputParams).maskParams(mask).listeners(Collections.singleton(listener)).create();
  }

  protected NgsParams getParams(final OutputStream out, final NgsMaskParams mask, final String subjects, final String queries, final int errorLimit, final boolean progress, final ListenerType listener, final boolean zip, final OutputFilter filter, final int topN) throws IOException {
    final File subjectsDir = FileHelper.createTempDirectory(mDir);
    final File queriesDir = FileHelper.createTempDirectory(mDir);
    final File hitsDir = FileHelper.createTempDirectory(mDir);
    //System.err.println("hitsDir=" + hitsDir);
    ReaderTestUtils.getReaderDNA(subjects, subjectsDir, null).close();
    ReaderTestUtils.getReaderDNA(queries, queriesDir, null).close();
    final SequenceParams subjectParams = SequenceParams.builder().directory(subjectsDir).mode(SequenceMode.UNIDIRECTIONAL).create();
    final SequenceParams queryParams = SequenceParams.builder().directory(queriesDir).create();
    final NgsFilterParams filterParams = NgsFilterParams.builder().outputFilter(filter).zip(zip).topN(topN).errorLimit(errorLimit).create();
    final File logDir = new File(hitsDir, "log");
    final NgsOutputParams outputParams = new NgsTestUtils.OverriddenNgsOutputParams(NgsTestUtils.OverriddenNgsOutputParams.builder().outStream(out).progress(progress).outputDir(logDir).filterParams(filterParams));
    return NgsParams.builder().stepSize(10).numberThreads(4).buildFirstParams(subjectParams).searchParams(queryParams).outputParams(outputParams).maskParams(mask).listeners(Collections.singleton(listener)).compressHashes(false).create();
  }

  @Override
  protected void check(final NgsMaskParams mask, final String subjects, final String queries, final String expected, final Long usageExp) throws Exception {
    check(mask, subjects, queries, expected, new String[] {
        "threads=4",
        "Performance statistics not available.",
    }, 4, usageExp);
  }

  @Override
  public void testC1Log() throws Exception {
    final ByteArrayOutputStream ba = new ByteArrayOutputStream();
    try (PrintStream pr = new PrintStream(ba)) {
      Diagnostic.setLogStream(pr);
      try {
        final ByteArrayOutputStream out = new ByteArrayOutputStream();
        final NgsParams params = getParams(out, new NgsMaskParamsGeneral(4, 0, 0, 1), SEQ_DNA_Y2, SEQ_DNA_Y1, 10/*errorlimit*/, true/*progress*/, ListenerType.NULL, false/*zip*/, OutputFilter.NONE, 1/*topN*/);
        try {
          execNgs(params);
        } finally {
          out.close();
        }
        assertEquals(NgsTestUtils.HEADER + EXPECTED_Y_Z, out.toString());
        //final File log = params.output().file("log");
        pr.flush();
        final String logs = ba.toString();
        //System.err.println(logs);
        //read counts
        assertTrue(logs.contains(" filter=NONE"));
        //marcher timers "
        TestUtils.containsAll(logs,
          " Thread Search parent Terminating ",
          " Thread Search 0 Start ",
          " Thread Search 0 Finish ",
          " Thread Search parent Finished ",
          " Start freeze job",
          " Finish freeze job"
        );
        //other timers
        //        assertTrue(logs.contains(" Timer Read_delay_subject "));
        if (params.numberThreads() > 1) {
          TestUtils.containsAll(logs,
            "Thread Search ", " Thread Search parent Terminating ", " Thread Search parent Finished ",
            " Thread Search 0 Scheduling", " Thread Search 0 Start ", " Thread Search 0 Finish ");
        }
        //deprecated timers to be removed
        TestUtils.containsAll(logs, " Timer Index_initialization ", " Timer Index_sort ", " Timer Index_pointer ", " Timer Index_position ", " Timer Index_bitVector ");
        //index statistics
        TestUtils.containsAll(logs, "Index[", "] statistics", "] search performance");
      } finally {
        Diagnostic.setLogStream();
      }
    } finally {
      ba.close();
    }
  }

  @Override
  public void testMultiCoreSplitting50() throws Exception {
    final String template = makeReads(">seq@" + LS + "ACGT" + LS, 51);
    final String reads = ""
        + ">read" + LS
        + "ACTG" + LS
        ;
    final String x = ""
        + "seq@" + TAB + "F" + TAB + "0" + TAB + "1" + TAB + "2" + TAB + "2" + LS
        + "seq@" + TAB + "R" + TAB + "0" + TAB + "1" + TAB + "2" + TAB + "2" + LS
        ;
    final String exp = makeReads(x, 51);

    final String[] logExp = {
        "threads=4",
        "Thread Search parent Finished",
        "Thread Search 0 Finish ",
    };
    check(new NgsMaskParamsGeneral(2, 1, 1, 1), reads, template, exp, logExp, 4, 4L);
  }

  @Override
  public void testMultiCoreSplitting40() throws Exception {
    final int max = 4;
    final String template = makeReads(">seq@" + LS + "ACGT" + LS, max);
    final String reads = ""
        + ">read" + LS
        + "ACTG" + LS
        ;
    final String x = ""
        + "seq@" + TAB + "F" + TAB + "0" + TAB + "1" + TAB + "2" + TAB + "2" + LS
        + "seq@" + TAB + "R" + TAB + "0" + TAB + "1" + TAB + "2" + TAB + "2" + LS
        ;
    final String exp = makeReads(x, max);

    final String[] logExp = {
        "threads=4",
        "Thread Search parent Finished",
        "Thread Search 0 Finish ",
    };
    check(new NgsMaskParamsGeneral(2, 1, 1, 1), reads, template, exp, logExp, 4, 4L);
  }

}
