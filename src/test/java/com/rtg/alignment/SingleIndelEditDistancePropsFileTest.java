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
