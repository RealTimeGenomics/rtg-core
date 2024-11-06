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
