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
public class CompleteGenomicsV2Machine extends CompleteGenomicsMachine {

  private static final int NUMBER_TRIES = 1000;
  private static final int READ_LENGTH = CgUtils.CG2_RAW_READ_LENGTH;
  private static final int MAX_DELETE = 20;  // Bring the right arm in from the end of the fragment to allow for some deletion.

  protected final double[] mOverlapDistribution2;

  /**
   * Constructs with seed and specific priors
   * @param params priors to use
   * @param randomSeed random seed
   */
  public CompleteGenomicsV2Machine(AbstractMachineErrorParams params, long randomSeed) {
    super(params, randomSeed);
    mOverlapDistribution2 = SimulationUtils.cumulativeDistribution(params.overlapDistribution2());
    if (mOverlapDistribution2.length != 8) {
      throw new IllegalArgumentException("Missing or incorrect distribution for CG V2 overlap");
    }
    mQualityBytes = new byte[READ_LENGTH];
    mReadBytes = new byte[READ_LENGTH];
    mWorkspace = new int[READ_LENGTH];
  }

  private int generateOverlapLength() {
    return SimulationUtils.chooseLength(mOverlapDistribution2, mSegmentRandom.nextDouble()) - 7;
  }

  @Override
  protected String generateRead(String id, int fragmentStart, byte[] data, int length, boolean forward, boolean leftArm) {
    final int startFrom;
    final int direction;
    final char frame;
    if (forward) {
      frame = 'F';
      direction = 1;
      if (leftArm) {
        startFrom = 0;
      } else {
        startFrom = length - (CgUtils.CG2_RAW_READ_LENGTH + MAX_DELETE);
      }
    } else {
      frame = 'R';
      direction = -1;
      if (leftArm) {
        startFrom = length - 1;
      } else {
        startFrom = CgUtils.CG2_RAW_READ_LENGTH + MAX_DELETE;
      }
    }
    for (int x = 0; x < NUMBER_TRIES; ++x) {
      resetCigar();
      int refPos = readBases(startFrom, data, leftArm ? length : startFrom + 1, direction, CgUtils.CG2_OVERLAP_POSITION);
      final int overlap = generateOverlapLength();
      refPos += overlap * direction;
      final int tLen = direction == 1 ? length - refPos : refPos;
      if (refPos < 0 || refPos >= length) {
        continue;
      }
      updateCigarWithPositiveOrNegativeSkip(overlap);
      refPos = readBases(refPos, data, tLen, direction, READ_LENGTH - CgUtils.CG2_OVERLAP_POSITION);
      //System.out.println("Generated overlap:" + overlap + " gap:" + gap);

      final int newStart = Math.min(startFrom, refPos - direction);
      if (!forward) {
        DNA.complementInPlace(mReadBytes, 0, mReadBytesUsed);
      }
      final String cigar = getCigar(!forward);
      //System.err.println(id + " fragmentStart: " + fragmentStart + " length: " + length + " forward: " + forward + " leftArm: " + leftArm + " startFrom:" + startFrom + " direction:" + direction);
      return formatReadName(id, frame, cigar, fragmentStart, newStart);
    }
    throw new NoTalkbackSlimException("Unable to generate a valid read with given priors in " + NUMBER_TRIES + " attempts");
  }

  @Override
  public MachineType machineType() {
    return MachineType.COMPLETE_GENOMICS_2;
  }
}
