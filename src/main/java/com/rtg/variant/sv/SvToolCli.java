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
package com.rtg.variant.sv;

import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;
import static com.rtg.util.cli.CommonFlagCategories.SENSITIVITY_TUNING;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.ParamsCli;
import com.rtg.launcher.SequenceParams;
import com.rtg.mode.SequenceMode;
import com.rtg.ngs.MapFlags;
import com.rtg.util.IORunnable;
import com.rtg.util.cli.CFlags;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.variant.sv.SvCliUtils.SvValidator;
import com.rtg.variant.sv.SvToolParams.SvToolParamsBuilder;

/**
 */
public class SvToolCli extends ParamsCli<SvToolParams> {

  private static final String BIN_SIZE = "Xbin-size";
  private static final String SV_SIMPLE = "simple-signals";
  private static final String FINE_STEP = "fine-step";
  private static final String CORRECTIONS_FLAG = "Xcorrections";
  private static final String HETEROZYGOUS_FLAG = "Xheterozygous";

  private static class SvToolValidator extends SvValidator {

    @Override
    public boolean isValid(final CFlags flags) {
      if (!super.isValid(flags)) {
        return false;
      }
      if ((Integer) flags.getValue(BIN_SIZE) < 1) {
        Diagnostic.error(ErrorType.EXPECTED_POSITIVE, BIN_SIZE);
        return false;
      }
      if ((Integer) flags.getValue(MapFlags.STEP_FLAG) < 1) {
        Diagnostic.error(ErrorType.EXPECTED_POSITIVE, MapFlags.STEP_FLAG);
        return false;
      }
      if ((Integer) flags.getValue(FINE_STEP) < 1) {
        Diagnostic.error(ErrorType.EXPECTED_POSITIVE, FINE_STEP);
        return false;
      }
      if ((Integer) flags.getValue(FINE_STEP) > (Integer) flags.getValue(MapFlags.STEP_FLAG)) {
        flags.setParseMessage("Parameter \"" + FINE_STEP + "\" should be smaller than or equal to parameter \"" + MapFlags.STEP_FLAG + "\"");
        return false;
      }
      return true;
    }
  }

  @Override
  public String moduleName() {
    return "sv";
  }

  @Override
  public String description() {
    return "find structural variant signals";
  }

  @Override
  protected void initFlags() {
    SvCliUtils.initCommonFlags(mFlags);
    SvCliUtils.initRelabelFlag(mFlags);
    mFlags.setDescription("Analyses SAM records to determine the location of structural variants.");
    mFlags.setValidator(new SvToolValidator());
    mFlags.registerOptional(SV_SIMPLE, "if set, also output simple signals").setCategory(INPUT_OUTPUT);
    mFlags.registerOptional('b', BIN_SIZE, Integer.class, CommonFlags.INT, "bin size used by simple signals", 10).setCategory(SENSITIVITY_TUNING);
    mFlags.registerOptional('s', MapFlags.STEP_FLAG, Integer.class, CommonFlags.INT, "step size", 100).setCategory(SENSITIVITY_TUNING);
    mFlags.registerOptional('f', FINE_STEP, Integer.class, CommonFlags.INT, "step size in interesting regions", 10).setCategory(SENSITIVITY_TUNING);
    //X flags
    mFlags.registerOptional(HETEROZYGOUS_FLAG, "if set, also include heterozygous bayesian models").setCategory(INPUT_OUTPUT);
    mFlags.registerOptional(CORRECTIONS_FLAG, File.class, CommonFlags.FILE, "file containing per position corrections").setCategory(INPUT_OUTPUT);
  }

  @Override
  protected SvToolParams makeParams() throws IOException {
    return makeParams(mFlags);
  }

  protected SvToolParams makeParams(CFlags flags) throws IOException {
    final SvToolParamsBuilder builder = SvToolParams.builder();
    SvCliUtils.populateCommonParams(builder, SequenceParams.builder().mode(SequenceMode.UNIDIRECTIONAL), flags);

    return builder.binSize((Integer) mFlags.getValue(BIN_SIZE))
        .stepSize((Integer) mFlags.getValue(MapFlags.STEP_FLAG))
        .fineStepSize((Integer) mFlags.getValue(FINE_STEP))
        .outputSimple(mFlags.isSet(SV_SIMPLE))
        .heterozygous(mFlags.isSet(HETEROZYGOUS_FLAG))
        .correctionsFile((File) mFlags.getValue(CORRECTIONS_FLAG))
        .create();
  }

  @Override
  protected File outputDirectory() {
    return (File) mFlags.getValue(CommonFlags.OUTPUT_FLAG);
  }

  @Override
  protected IORunnable task(final SvToolParams params, final OutputStream out) throws IOException {
    return new SvToolTask(params, out);
  }

}
