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

import java.io.File;
import java.io.IOException;

import com.rtg.launcher.AbstractNanoTest;
import com.rtg.mode.ProteinScoringMatrix;
import com.rtg.ngs.NgsParams;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;

/**
 * Test class.
 */
public class TopNProteinOutputProcessorTest extends AbstractNanoTest {

  public void testTopN() throws IOException, InvalidParamsException {
    try (final TestDirectory tmp = new TestDirectory("proteintopn")) {
      final NgsParams params = TopEqualProteinOutputProcessorTest.createParams(tmp, 3);
      try (final TopNProteinOutputProcessor p = new TopNProteinOutputProcessor(params, null)) {
        p.writeResult(getProteinRes(0, -42, params));
        p.writeResult(getProteinRes(1, 42, params));
        p.writeResult(getProteinRes(1, -42, params));
        p.writeResult(getProteinRes(2, -42, params));
        p.writeResult(getProteinRes(2, 42, params));
        p.writeResult(getProteinRes(3, -42, params));
        p.writeResult(getProteinRes(3, -42, params));
        p.writeResult(getProteinRes(4, -42, params));
        p.writeResult(getProteinRes(4, -42, params));
        p.writeResult(getProteinRes(4, -42, params));
        p.writeResult(getProteinRes(5, -42, params));
        p.writeResult(getProteinRes(5, -42, params));
        p.writeResult(getProteinRes(5, -42, params));
        p.writeResult(getProteinRes(5, -42, params));
        p.writeResult(getProteinRes(6, 9, params));
        p.writeResult(getProteinRes(6, 2, params));
        p.writeResult(getProteinRes(6, 2, params));
        p.writeResult(getProteinRes(6, 3, params));
        p.writeResult(getProteinRes(6, 3, params));
        p.finish();
      }
      mNano.check("test-topn.txt", TestUtils.sanitizeTsvHeader(StringUtils.grepMinusV(FileUtils.fileToString(new File(tmp, "alignments.tsv")), "^#.*-ID")));
    }
  }

  public void testTopN250() throws IOException, InvalidParamsException {
    try (final TestDirectory tmp = new TestDirectory("proteintopn")) {
      final NgsParams params = TopEqualProteinOutputProcessorTest.createParams(tmp, 250);
      try (final TopNProteinOutputProcessor p = new TopNProteinOutputProcessor(params, null)) {
        for (int k = 0; k < 250; ++k) {
          p.writeResult(getProteinRes(0, -47, params));
        }
        for (int k = 0; k < 251; ++k) {
          p.writeResult(getProteinRes(1, 47, params));
        }
        p.finish();
      }
      mNano.check("test-topn250.txt", TestUtils.sanitizeTsvHeader(StringUtils.grepMinusV(FileUtils.fileToString(new File(tmp, "alignments.tsv")), "^#.*-ID")));
    }
  }

  private ProteinAlignmentResult getProteinRes(int readid, final int score, NgsParams params) throws InvalidParamsException, IOException {
    final SharedProteinResources resx = new SharedProteinResources(new ProteinScoringMatrix(), params.searchParams().reader(), params.buildFirstParams().reader(), false);
    return new ProteinAlignmentResult(resx, 0, readid * 6, null, 0, true) {
        @Override
        int alignmentScore() {
          return score;
        }
      };
  }

}
