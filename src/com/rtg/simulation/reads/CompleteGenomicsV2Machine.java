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
import com.rtg.util.PortableRandom;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.variant.AbstractMachineErrorParams;

/**
 * Generator for Complete Genomics paired end reads.
 *
 */
public class CompleteGenomicsV2Machine extends CompleteGenomicsMachine {

  private static final int NUMBER_TRIES = 1000;
  private static final int READ_LENGTH = CgUtils.CG2_RAW_LENGTH;

  protected final PortableRandom mFrameRandom;
  protected final PortableRandom mSegmentRandom;
  protected final double[] mOverlapDistribution2;

  /**
   * Constructs with seed and specific priors
   * @param params priors to use
   * @param randomSeed random seed
   */
  public CompleteGenomicsV2Machine(AbstractMachineErrorParams params, long randomSeed) {
    super(params, randomSeed);
    mOverlapDistribution2 = SimulationUtils.cumulativeDistribution(params.overlapDistribution2());
    if (mOverlapDistribution2 == null || mOverlapDistribution2.length != 8) {
      throw new IllegalArgumentException("Missing or incorrect distribution for CG V2 overlap");
    }
    mFrameRandom = new PortableRandom(randomSeed);
    mSegmentRandom = new PortableRandom(randomSeed * 3);
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
    } else {
      frame = 'R';
    }
    if (forward) {
      direction = 1;
      startFrom = 0;
    } else {
      direction = -1;
      startFrom = length - 1;
    }
    for (int x = 0; x < NUMBER_TRIES; x++) {
      resetCigar();
      int refPos = readBases(startFrom, data, length, direction, 10);
      final int overlap = generateOverlapLength();
      refPos += overlap * direction;
      final int tLen = direction == 1 ? length - refPos : length - startFrom + refPos;
      if (refPos < 0 || refPos >= length) {
        continue;
      }
      updateCigarWithPositiveOrNegativeSkip(overlap);
      refPos = readBases(refPos, data, tLen, direction, 19);
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
}
