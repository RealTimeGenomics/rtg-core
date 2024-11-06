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
package com.rtg.alignment;

import static com.rtg.util.StringUtils.LS;

import java.io.File;
import java.io.IOException;

import com.rtg.ngs.NgsParams;
import com.rtg.ngs.NgsParamsBuilder;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;

/**
 * Load a props file that is the equivalent of the affine penalties.
 */
public class SingleIndelEditDistancePropsFileTest extends SingleIndelEditDistanceTest {
  private static final String PROPS_TEMPLATE = ""
                                         + "error_snp_penalty=%d" + LS
                                         + "error_unknowns_penalty=%d" + LS
                                         + "error_del_penalty=%d" + LS
                                         + "error_ins_penalty=%d" + LS
                                         + "error_ins_penalty_extension_slope=%f" + LS
                                         + "error_del_penalty_extension_slope=%f" + LS
      ;
  @Override
  protected UnidirectionalEditDistance getEditDistanceInstance(int gapOpenPenalty, int gapExtendPenalty, int substitutionPenalty, int unknownsPenalty) {
    try {
      try (TestDirectory dir = new TestDirectory()) {
        final File testProperties = createPropsFile(gapOpenPenalty, gapExtendPenalty, substitutionPenalty, unknownsPenalty, dir);
        final NgsParams params = new NgsParamsBuilder().singleIndelPenalties(testProperties.getPath()).create();
        return new SingleIndelEditDistance(params, 1000);
      }
    } catch (IOException e) {
      throw new RuntimeException("Can't create test directory");
    }
  }

  private File createPropsFile(int gapOpenPenalty, int gapExtendPenalty, int substitutionPenalty, int unknownsPenalty, TestDirectory dir) throws IOException {
    final File testProperties = new File(dir, "props");
    FileUtils.stringToFile(String.format(PROPS_TEMPLATE
        , substitutionPenalty
        , unknownsPenalty
        , gapOpenPenalty + gapExtendPenalty
        , gapOpenPenalty + gapExtendPenalty
        , (double) gapExtendPenalty
        , (double) gapExtendPenalty
    )
        , testProperties);
    return testProperties;
  }
}
