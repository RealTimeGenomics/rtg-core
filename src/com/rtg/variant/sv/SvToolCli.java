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
import com.rtg.util.IORunnable;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.cli.CFlags;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.variant.sv.SvCliUtils.SvValidator;
import com.rtg.variant.sv.SvToolParams.SvToolParamsBuilder;

/**
 */
public class SvToolCli extends ParamsCli<SvToolParams> {

  private static final String MODULE_NAME = "sv";
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
      if ((Integer) flags.getValue(CommonFlags.STEP_FLAG) < 1) {
        Diagnostic.error(ErrorType.EXPECTED_POSITIVE, CommonFlags.STEP_FLAG);
        return false;
      }
      if ((Integer) flags.getValue(FINE_STEP) < 1) {
        Diagnostic.error(ErrorType.EXPECTED_POSITIVE, FINE_STEP);
        return false;
      }
      if ((Integer) flags.getValue(FINE_STEP) > (Integer) flags.getValue(CommonFlags.STEP_FLAG)) {
        flags.error("Parameter \"" + FINE_STEP + "\" should be smaller than or equal to parameter \"" + CommonFlags.STEP_FLAG + "\"");
        return false;
      }
      return true;
    }
  }

  @Override
  public String moduleName() {
    return MODULE_NAME;
  }

  @Override
  protected void initFlags() {
    initFlags(mFlags);
  }

  protected static void initFlags(CFlags flags) {
    SvCliUtils.initCommonFlags(flags);
    SvCliUtils.initRelabelFlag(flags);
    flags.setDescription("Analyses SAM records to determine the location of structural variants.");
    flags.setValidator(new SvToolValidator());
    flags.registerOptional(SV_SIMPLE, "if set, also output simple signals").setCategory(INPUT_OUTPUT);
    flags.registerOptional('b', BIN_SIZE, Integer.class, "INT", "bin size used by simple signals", 10).setCategory(SENSITIVITY_TUNING);
    flags.registerOptional('s', CommonFlags.STEP_FLAG, Integer.class, "INT", "step size", 100).setCategory(SENSITIVITY_TUNING);
    flags.registerOptional('f', FINE_STEP, Integer.class, "INT", "step size in interesting regions", 10).setCategory(SENSITIVITY_TUNING);
    //X flags
    flags.registerOptional(HETEROZYGOUS_FLAG, "if set, also include heterozygous bayesian models").setCategory(INPUT_OUTPUT);
    flags.registerOptional(CORRECTIONS_FLAG, File.class, "FILE", "file containing per position corrections").setCategory(INPUT_OUTPUT);
  }

  @Override
  protected SvToolParams makeParams() throws InvalidParamsException, IOException {
    return makeParams(mFlags);
  }

  protected SvToolParams makeParams(CFlags flags) throws InvalidParamsException, IOException {
    final SvToolParamsBuilder builder = SvToolParams.builder();
    SvCliUtils.populateCommonParams(builder, SequenceParams.builder().mode(SequenceMode.UNIDIRECTIONAL), flags);

    return builder.binSize((Integer) mFlags.getValue(BIN_SIZE))
        .stepSize((Integer) mFlags.getValue(CommonFlags.STEP_FLAG))
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
