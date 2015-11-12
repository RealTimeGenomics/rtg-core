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

import java.io.IOException;

import com.rtg.util.InvalidParamsException;
import com.rtg.util.machine.MachineType;
import com.rtg.variant.AbstractMachineErrorParams;
import com.rtg.variant.MachineErrorParamsBuilder;

/**
 * 454 single end machine
 */
public class FourFiveFourSingleEndMachine extends SingleEndRandomLengthMachine {

  /**
   * Construct with given priors and seed
   * @param params priors
   * @param randomSeed seed for random number generation
   */
  public FourFiveFourSingleEndMachine(AbstractMachineErrorParams params, long randomSeed) {
    super(params, randomSeed);
  }

  /**
   * Construct with given random seed and default 454 priors
   * @param randomSeed seed for random number generation
   * @throws InvalidParamsException If priors fail to load
   * @throws IOException whenever
   */
  public FourFiveFourSingleEndMachine(long randomSeed) throws InvalidParamsException, IOException {
    this(new MachineErrorParamsBuilder().errors("ls454_se").create(), randomSeed);
  }

  @Override
  public MachineType machineType() {
    return MachineType.FOURFIVEFOUR_SE;
  }
}
