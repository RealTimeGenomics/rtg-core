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

import com.rtg.reader.PrereadType;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.PortableRandom;
import com.rtg.variant.AbstractMachineErrorParams;
import com.rtg.variant.MachineErrorParamsBuilder;

/**
 * super class for both single and paired end Illumina machines
 */
public abstract class AbstractIlluminaMachine extends AbstractMachine {

  protected final PortableRandom mFrameRandom;


  /**
   * Constructs with seed and specific priors
   * @param params priors to use
   * @param randomSeed random seed
   */
  public AbstractIlluminaMachine(AbstractMachineErrorParams params, long randomSeed) {
    super(params);
    mFrameRandom = new PortableRandom(randomSeed);
  }

  /**
   * Constructs with seed and default Illumina priors
   * @param randomSeed random seed
   * @throws InvalidParamsException if fails to construct priors
   * @throws IOException whenever
   */
  public AbstractIlluminaMachine(long randomSeed) throws InvalidParamsException, IOException {
    this(new MachineErrorParamsBuilder().errors("illumina").create(), randomSeed);
  }

  protected String generateRead(String id, int fragmentStart, byte[] data, int length, boolean forward, int readLength) {
    final int startFrom;
    final int direction;
    if (forward) {
      startFrom = 0;
      direction = 1;
    } else {
      startFrom = length - 1;
      direction = -1;
    }
    final int localStart = process(startFrom, data, length, direction, readLength);
    final String cigar = getCigar(!forward, localStart, length, readLength);
    return formatReadName(id, forward ? 'F' : 'R', cigar, fragmentStart, localStart);
  }

  @Override
  public abstract void processFragment(String id, int fragmentStart, byte[] data, int length) throws IOException;


  @Override
  public abstract boolean isPaired();

  @Override
  public PrereadType machineType() {
    return PrereadType.SOLEXA;
  }

}
