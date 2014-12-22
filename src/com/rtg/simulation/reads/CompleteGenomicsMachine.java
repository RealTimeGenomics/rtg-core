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

import com.rtg.alignment.ActionsHelper;
import com.rtg.mode.DNA;
import com.rtg.reader.PrereadType;
import com.rtg.simulation.SimulationUtils;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.PortableRandom;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.variant.AbstractMachineErrorParams;
import com.rtg.variant.MachineErrorParamsBuilder;

/**
 * Generator for Complete Genomics paired end reads.
 *
 */
public class CompleteGenomicsMachine extends AbstractMachine {

  private static final int NUMBER_TRIES = 1000;
  private static final int READ_LENGTH = 35;

  protected final PortableRandom mFrameRandom;
  protected final PortableRandom mSegmentRandom;
  protected final double[] mOverlapDistribution;
  protected final double[] mGapDistribution;
  protected final double[] mSmallGapDistribution;

  /**
   * Constructs with seed and specific priors
   * @param params priors to use
   * @param randomSeed random seed
   */
  public CompleteGenomicsMachine(AbstractMachineErrorParams params, long randomSeed) {
    super(params);
    mOverlapDistribution = accumDistribution(params.overlapDistribution());
    mGapDistribution = accumDistribution(params.gapDistribution());
    mSmallGapDistribution = accumDistribution(params.smallGapDistribution());
    mFrameRandom = new PortableRandom(randomSeed);
    mSegmentRandom = new PortableRandom(randomSeed * 3);
    setBuffers();
  }

  /**
   * Constructs with seed and default Illumina priors
   * @param randomSeed random seed
   * @throws InvalidParamsException if fails to construct priors
   * @throws IOException whenever
   */
  public CompleteGenomicsMachine(long randomSeed) throws InvalidParamsException, IOException {
    this(new MachineErrorParamsBuilder().errors("complete").create(), randomSeed);
  }


  private static double[] accumDistribution(final double[] dist) {
    final double[] accumDist = new double[dist.length];
    double accum = 0.0;
    for (int i = 0; i < dist.length; i++) {
      accum += dist[i];
      accumDist[i] = accum;
    }
    return accumDist;
  }

  private void setBuffers() {
    mQualityBytes = new byte[READ_LENGTH];
    mReadBytes = new byte[READ_LENGTH];
    mWorkspace = new int[READ_LENGTH];
  }

  private void reverse() {
    for (int left = 0, right = mReadBytesUsed - 1; left < right; left++, right--) {
      // exchange the first and last
      final byte temp = mReadBytes[left];
      mReadBytes[left] = mReadBytes[right];
      mReadBytes[right] = temp;
    }
  }

  private int generateOverlapLength() {
    if ((mOverlapDistribution == null) || (mOverlapDistribution.length != 5)) {
      return -2;
    }
    return SimulationUtils.chooseLength(mOverlapDistribution, mSegmentRandom.nextDouble()) - 4;
  }

  private int generateGapLength() {
    if ((mGapDistribution == null) || (mGapDistribution.length != 5)) {
      return 6;
    }
    return SimulationUtils.chooseLength(mGapDistribution, mSegmentRandom.nextDouble()) + 4;
  }

  private int generateSmallGapLength() {
    if ((mSmallGapDistribution == null) || (mSmallGapDistribution.length != 4)) {
      System.err.println("wrong");
      return 0;
    }
    return SimulationUtils.chooseLength(mSmallGapDistribution, mSegmentRandom.nextDouble());
  }

  private void updateCigarWithPositiveOrNegativeSkip(final int skip) {
    if (skip > 0) {
      // N operation
      addCigarState(skip, ActionsHelper.CG_GAP_IN_READ);
    } else if (skip < 0) {
      // B operation
      addCigarState(-skip, ActionsHelper.CG_OVERLAP_IN_READ);
    }
  }

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
    Diagnostic.developerLog(id + " fragmentStart: " + fragmentStart + " length: " + length + " forward: " + forward + " leftArm: " + leftArm);
    throw new NoTalkbackSlimException("Unable to generate a valid read with given priors in " + NUMBER_TRIES + " attempts");
  }

  @Override
  public void processFragment(String id, int fragmentStart, byte[] data, int length) throws IOException {
    reseedErrorRandom(mFrameRandom.nextLong());
    final boolean forwardFrame = mFrameRandom.nextBoolean();

    final String nameLeft = generateRead(id, fragmentStart, data, length, forwardFrame, true);

    mReadWriter.writeLeftRead(nameLeft, mReadBytes, mQualityBytes, READ_LENGTH);
    mResidueCount += READ_LENGTH;

    final String nameRight = generateRead(id, fragmentStart, data, length, forwardFrame, false);
    mReadWriter.writeRightRead(nameRight, mReadBytes, mQualityBytes, READ_LENGTH);
    mResidueCount += READ_LENGTH;
  }


  @Override
  public boolean isPaired() {
    return true;
  }

  @Override
  public PrereadType machineType() {
    return PrereadType.CG;
  }

}
