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
import com.rtg.util.PortableRandom;
import com.rtg.variant.AbstractMachineErrorParams;

/**
 * 454 paired end machine
 */
public class FourFiveFourPairedEndMachine extends AbstractMachine {

  private static final int NUMBER_TRIES = 1000;
  private static final int MIN_SIDE_LENGTH = 1; // Ensure the smallest side has at least this many bases

  private int mMinPairSize;
  private int mMaxPairSize;
  protected final PortableRandom mPairSizeRandom;
  protected final PortableRandom mPairPositionRandom;
  protected final PortableRandom mFrameRandom;

  /**
   * Construct with given priors and seed
   * @param params priors
   * @param randomSeed seed for random number generation
   */
  public FourFiveFourPairedEndMachine(AbstractMachineErrorParams params, long randomSeed) {
    super(params);
    mPairSizeRandom = new PortableRandom(randomSeed);
    mPairPositionRandom = new PortableRandom();
    mFrameRandom = new PortableRandom();
  }

  protected void reseedErrorRandom() {
    final long seed = mPairSizeRandom.nextLong();
    mPairPositionRandom.setSeed(seed + 1);
    mFrameRandom.setSeed(seed + 2);
    super.reseedErrorRandom(seed + 2);
  }

  /**
   * Set the minimum size of the total lengths of both sides of a read pair
   * @param size the size
   */
  public void setMinPairSize(int size) {
    mMinPairSize = size;
  }

  /**
   * Set the maximum size of the total length of a read pair
   * @param size the size
   */
  public void setMaxPairSize(int size) {
    mMaxPairSize = size;
  }

  @Override
  public boolean isPaired() {
    return true;
  }

  @Override
  public PrereadType machineType() {
    return PrereadType.UNKNOWN;
  }

  void updateWorkingSpace(int length) {
    if (mReadBytes.length < length) {
      mReadBytes = new byte[length];
      mWorkspace = new int[length];
      mQualityBytes = new byte[length];
    }
  }

  @Override
  public void processFragment(String id, int fragmentStart, byte[] data, int length) throws IOException {
    updateWorkingSpace(length);
    reseedErrorRandom();
    final double mid = (mMaxPairSize + mMinPairSize) * 0.5 + 0.5;
    final double width = (mMaxPairSize - mMinPairSize) * 0.25; // 2 std devs per side
    int pairLength = 0;
    for (int x = 0; x < NUMBER_TRIES; x++) {
      pairLength = (int) (mPairSizeRandom.nextGaussian() * width + mid);
      if ((pairLength >= mMinPairSize) && (pairLength <= mMaxPairSize)) {
        break;
      }
    }

    //this is the position on the subFragment that we cross the adapter
    final int pairPosition = mPairPositionRandom.nextInt(pairLength - (2 * MIN_SIDE_LENGTH)) + MIN_SIDE_LENGTH;

    final boolean forward = mFrameRandom.nextBoolean();
    int pos;
    if (forward) {
      pos = process(length - pairLength + pairPosition, data, pairLength - pairPosition, 1, length);
    } else {
      pos = process(pairPosition, data, pairPosition, -1, length);
    }
    String cigar = getCigar(!forward, pos, length, mReadBytesUsed);
    String name = formatReadName(id, forward ? 'F' : 'R', cigar, fragmentStart, pos);
    mReadWriter.writeLeftRead(name, mReadBytes, mQualityBytes, mReadBytesUsed);
    mResidueCount += mReadBytesUsed;
      //System.out.println("length=" + length + " pairLength: " + pairLength + " pp: " + pairPosition + " pos=" + pos + " fs=" + fragmentStart + " name" + name);

    if (forward) {
      pos = process(0, data, pairPosition, 1, length);
    } else {
      pos = process(length - 1, data, pairLength - pairPosition, -1, length);
    }
    cigar = getCigar(!forward, pos, length, mReadBytesUsed);
    name = formatReadName(id, forward ? 'F' : 'R', cigar, fragmentStart, pos);
    mReadWriter.writeRightRead(name, mReadBytes, mQualityBytes, mReadBytesUsed);
    mResidueCount += mReadBytesUsed;
  }
}
