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
