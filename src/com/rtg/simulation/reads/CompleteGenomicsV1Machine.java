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

import com.rtg.mode.DNA;
import com.rtg.reader.CgUtils;
import com.rtg.simulation.SimulationUtils;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.machine.MachineType;
import com.rtg.variant.AbstractMachineErrorParams;

/**
 * Generator for Complete Genomics paired end reads.
 *
 */
public class CompleteGenomicsV1Machine extends CompleteGenomicsMachine {

  private static final int NUMBER_TRIES = 1000;
  private static final int READ_LENGTH = CgUtils.CG_RAW_READ_LENGTH;

  protected final double[] mOverlapDistribution;
  protected final double[] mGapDistribution;
  protected final double[] mSmallGapDistribution;

  /**
   * Constructs with seed and specific priors
   * @param params priors to use
   * @param randomSeed random seed
   */
  public CompleteGenomicsV1Machine(AbstractMachineErrorParams params, long randomSeed) {
    super(params, randomSeed);
    mOverlapDistribution = SimulationUtils.cumulativeDistribution(params.overlapDistribution());
    mGapDistribution = SimulationUtils.cumulativeDistribution(params.gapDistribution());
    mSmallGapDistribution = SimulationUtils.cumulativeDistribution(params.smallGapDistribution());
    if (mOverlapDistribution.length != 5) {
      throw new IllegalArgumentException("Missing or incorrect distribution for CG V1 overlap");
    }
    if (mGapDistribution.length != 5) {
      throw new IllegalArgumentException("Missing or incorrect distribution for CG V1 gap");
    }
    if (mSmallGapDistribution.length != 4) {
      throw new IllegalArgumentException("Missing or incorrect distribution for CG V1 small gap");
    }
    mQualityBytes = new byte[READ_LENGTH];
    mReadBytes = new byte[READ_LENGTH];
    mWorkspace = new int[READ_LENGTH];
  }

  private int generateOverlapLength() {
    return SimulationUtils.chooseLength(mOverlapDistribution, mSegmentRandom.nextDouble()) - 4;
  }

  private int generateGapLength() {
    return SimulationUtils.chooseLength(mGapDistribution, mSegmentRandom.nextDouble()) + 4;
  }

  private int generateSmallGapLength() {
    return SimulationUtils.chooseLength(mSmallGapDistribution, mSegmentRandom.nextDouble());
  }

  @Override
  protected String generateRead(String id, int fragmentStart, byte[] data, int length, boolean forward, boolean leftArm) {
    final int startFrom;
    final int direction;
    final char frame;
    if (forward) {
      frame = 'F';
    } else {
      frame = 'R';
    }
    if (leftArm ^ !forward) {
      direction = 1;
      startFrom = 0;
    } else {
      direction = -1;
      startFrom = length - 1;
    }
    for (int x = 0; x < NUMBER_TRIES; x++) {
      resetCigar();
      int refPos = readBases(startFrom, data, length, direction, 5);
      final int overlap = generateOverlapLength();
      refPos += overlap * direction;
      int tLen = direction == 1 ? length - refPos : length - startFrom + refPos;
      if (refPos < 0 || refPos >= length) {
        continue;
      }
      updateCigarWithPositiveOrNegativeSkip(overlap);
      refPos = readBases(refPos, data, tLen, direction, 10);

      // Currently hardcoded gap of 0
      final int smallgap = generateSmallGapLength();
      refPos += smallgap * direction;
      tLen = direction == 1 ? length - refPos : length - startFrom + refPos;
      if (refPos < 0 || refPos >= length) {
        continue;
      }
      updateCigarWithPositiveOrNegativeSkip(smallgap);
      refPos = readBases(refPos, data, tLen, direction, 10);

      final int gap = generateGapLength();
      refPos += gap * direction;
      tLen = direction == 1 ? length - refPos : length - startFrom + refPos;
      if (refPos < 0 || refPos >= length) {
        continue;
      }
      updateCigarWithPositiveOrNegativeSkip(gap);
      refPos = readBases(refPos, data, tLen, direction, 10);
      //System.out.println("Generated overlap:" + overlap + " gap:" + gap);

      final int newStart = Math.min(startFrom, refPos - direction);
      if (forward ^ direction == 1) {
        reverse();
      }
      if (!forward) {
        DNA.complementInPlace(mReadBytes, 0, mReadBytesUsed);
      }
      final String cigar = getCigar(direction == -1, newStart, length, READ_LENGTH);
      return formatReadName(id, frame, cigar, fragmentStart, newStart);
    }
    //Diagnostic.developerLog(id + " fragmentStart: " + fragmentStart + " length: " + length + " forward: " + forward + " leftArm: " + leftArm);
    throw new NoTalkbackSlimException("Unable to generate a valid read with given priors in " + NUMBER_TRIES + " attempts");
  }

  @Override
  public MachineType machineType() {
    return MachineType.COMPLETE_GENOMICS;
  }
}
