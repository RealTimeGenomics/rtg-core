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

import org.junit.runner.RunWith;
import org.junit.runners.Suite;

import com.rtg.launcher.ModuleParamsTest;
import com.rtg.launcher.OutputModuleParamsTest;
import com.rtg.metagenomics.SpeciesParamsTest;
import com.rtg.ngs.NgsParamsTest;
import com.rtg.sam.MappedParamsTest;
import com.rtg.sam.SingleMappedParamsTest;
import com.rtg.variant.VariantParamsTest;
import com.rtg.variant.cnv.CnvProductParamsTest;
import com.rtg.variant.coverage.CoverageParamsTest;
import com.rtg.variant.sv.SvParamsTest;
import com.rtg.variant.sv.SvToolParamsTest;
import com.rtg.variant.sv.discord.DiscordantToolParamsTest;
import com.rtg.vcf.eval.VcfEvalParamsTest;

/**
 */
@RunWith(Suite.class)
@Suite.SuiteClasses({
  ModuleParamsTest.class,
  OutputModuleParamsTest.class,
  SpeciesParamsTest.class,
  NgsParamsTest.class,
  MappedParamsTest.class,
  SingleMappedParamsTest.class,
  CoverageParamsTest.class,
  VariantParamsTest.class,
  CnvProductParamsTest.class,
  VcfEvalParamsTest.class,
  SvParamsTest.class,
  SvToolParamsTest.class,
  DiscordantToolParamsTest.class,
})
public class AllTestParamsSuite {
  // required empty suite class
}
