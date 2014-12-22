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
package com.rtg.simulation.reads;

import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;
import static com.rtg.util.cli.CommonFlagCategories.UTILITY;

import com.rtg.util.cli.CFlags;
import com.rtg.util.machine.MachineType;

/**
 * Thin CLI layer for Complete genomics simulator
 */
public class CgSimCli extends ReadSimCli {

  static final String XMACHINE_ERROR_PRIORS = "Xmachine-errors";

  protected static class CgSimValidator extends ReadSimValidator {
    @Override
    protected boolean checkMachines(CFlags cflags) {
      return true;
    }
  }

  @Override
  public String moduleName() {
    return "cgsim";
  }

  @Override
  protected String getPriorsNameFlagValue() {
    if (mFlags.isSet(XMACHINE_ERROR_PRIORS)) {
      return (String) mFlags.getValue(XMACHINE_ERROR_PRIORS);
    }
    return null;
  }

  @Override
  protected MachineType getMachineType() {
    return MachineType.COMPLETE_GENOMICS;
  }

  @Override
  protected void initFlags() {
    super.initFlags();
    mFlags.setDescription("Simulate mutations in Complete Genomics reads.");
    mFlags.setCategories(UTILITY, new String[]{INPUT_OUTPUT, CAT_FRAGMENTS, UTILITY});
    mFlags.setValidator(new CgSimValidator());

  }
  @Override
  protected void initMachineFlags() {
    mFlags.registerOptional('M', MAX_FRAGMENT, Integer.class, "int", "maximum fragment size", 500).setCategory(CAT_FRAGMENTS);
    mFlags.registerOptional('m', MIN_FRAGMENT, Integer.class, "int", "minimum fragment size", 350).setCategory(CAT_FRAGMENTS);
    mFlags.registerOptional('E', XMACHINE_ERROR_PRIORS, String.class, "string", "override default machine error priors").setCategory(UTILITY);
  }
}
