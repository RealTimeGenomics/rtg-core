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
package com.rtg.util.test.params;

import com.rtg.launcher.ModuleParamsTest;
import com.rtg.launcher.OutputModuleParamsTest;
import com.rtg.metagenomics.SpeciesParamsTest;
import com.rtg.ngs.NgsParamsTest;
import com.rtg.sam.MappedParamsTest;
import com.rtg.sam.SingleMappedParamsTest;
import com.rtg.variant.VariantParamsTest;
import com.rtg.variant.cnv.CnvProductParamsTest;
import com.rtg.variant.coverage.CoverageParamsTest;
import com.rtg.variant.eval.VcfEvalParamsTest;
import com.rtg.variant.sv.SvParamsTest;
import com.rtg.variant.sv.SvToolParamsTest;
import com.rtg.variant.sv.discord.DiscordantToolParamsTest;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 */
public class AllTestParams extends TestSuite {

  /**
   * @return a test suite for all test classes that use <code>TestParams</code>
   */
  public static Test suite() {
    final TestSuite suite = new TestSuite("com.rtg.util.test.params");
    suite.addTestSuite(ModuleParamsTest.class);
    suite.addTestSuite(OutputModuleParamsTest.class);
    suite.addTestSuite(SpeciesParamsTest.class);
    suite.addTestSuite(NgsParamsTest.class);
    suite.addTestSuite(MappedParamsTest.class);
    suite.addTestSuite(SingleMappedParamsTest.class);
    suite.addTestSuite(CoverageParamsTest.class);
    suite.addTestSuite(VariantParamsTest.class);
    suite.addTestSuite(CnvProductParamsTest.class);
    suite.addTestSuite(VcfEvalParamsTest.class);
    suite.addTestSuite(SvParamsTest.class);
    suite.addTestSuite(SvToolParamsTest.class);
    suite.addTestSuite(DiscordantToolParamsTest.class);
    return suite;
  }

}
