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

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.reader.PrereadType;
import com.rtg.util.PortableRandom;
import com.rtg.variant.AbstractMachineErrorParams;

/**
 * single end machine with random read lengths
 */
@TestClass(value = {"com.rtg.simulation.reads.FourFiveFourSingleEndMachineTest", "com.rtg.simulation.reads.IonTorrentSingleEndMachineTest"})
public class SingleEndRandomLengthMachine extends AbstractMachine {

  private static final int NUMBER_TRIES = 1000;

  private int mMinSize;
  private int mMaxSize;
  protected final PortableRandom mReadSizeRandom;
//  protected final PortableRandom mPairPositionRandom;
  protected final PortableRandom mFrameRandom;

  /**
   * Construct with given priors and seed
   * @param params priors
   * @param randomSeed seed for random number generation
   */
  public SingleEndRandomLengthMachine(AbstractMachineErrorParams params, long randomSeed) {
    super(params);
    mReadSizeRandom = new PortableRandom(randomSeed);
//    mPairPositionRandom = new PortableRandom();
    mFrameRandom = new PortableRandom();
  }

  protected void reseedErrorRandom() {
    final long seed = mReadSizeRandom.nextLong();
//    mPairPositionRandom.setSeed(seed + 1);
    mFrameRandom.setSeed(seed + 2);
    super.reseedErrorRandom(seed + 2);
  }

  /**
   * Set the minimum size of the total lengths of both sides of a read
   * @param size the size
   */
  public void setMinSize(int size) {
    mMinSize = size;
  }

  /**
   * Set the maximum size of the total length of a read
   * @param size the size
   */
  public void setMaxSize(int size) {
    mMaxSize = size;
  }

  @Override
  public boolean isPaired() {
    return false;
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
    final double mid = (mMaxSize + mMinSize) * 0.5 + 0.5;
    final double width = (mMaxSize - mMinSize) * 0.25; // 2 std devs per side
    int readLength = 0;
    for (int x = 0; x < NUMBER_TRIES; x++) {
      readLength = (int) (mReadSizeRandom.nextGaussian() * width + mid);
      if ((readLength >= mMinSize) && (readLength <= mMaxSize)) {
        break;
      }
    }
    final boolean forward = mFrameRandom.nextBoolean();
    final int pos;
    if (forward) {
      pos = processBackwards(length - 1, data, length, -1, readLength);
    } else {
      pos = processBackwards(0, data, length, 1, readLength);
    }
    final String cigar = getCigar(forward, pos, length, mReadBytesUsed);  //forward rather than !forward because we used processBackwards
    final String name = formatReadName(id, forward ? 'F' : 'R', cigar, fragmentStart, pos);
    mReadWriter.writeRead(name, mReadBytes, mQualityBytes, mReadBytesUsed);
    mResidueCount += mReadBytesUsed;
  }
}
