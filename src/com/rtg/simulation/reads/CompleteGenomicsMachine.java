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
import com.rtg.alignment.ActionsHelper;
import com.rtg.reader.PrereadType;
import com.rtg.util.PortableRandom;
import com.rtg.variant.AbstractMachineErrorParams;

/**
 * Generator for Complete Genomics paired end reads.
 *
 */
@TestClass("com.rtg.simulation.reads.CompleteGenomicsV1MachineTest")
public abstract class CompleteGenomicsMachine extends AbstractMachine {

  protected final PortableRandom mFrameRandom;
  protected final PortableRandom mSegmentRandom;

  /**
   * Constructs with seed and specific priors
   * @param params priors to use
   * @param randomSeed random seed
   */
  public CompleteGenomicsMachine(AbstractMachineErrorParams params, long randomSeed) {
    super(params);
    mFrameRandom = new PortableRandom(randomSeed);
    mSegmentRandom = new PortableRandom(randomSeed * 3);
  }

  protected void reverse() {
    for (int left = 0, right = mReadBytesUsed - 1; left < right; left++, right--) {
      // exchange the first and last
      final byte temp = mReadBytes[left];
      mReadBytes[left] = mReadBytes[right];
      mReadBytes[right] = temp;
    }
  }

  protected void updateCigarWithPositiveOrNegativeSkip(final int skip) {
    if (skip > 0) {
      // N operation
      addCigarState(skip, ActionsHelper.CG_GAP_IN_READ);
    } else if (skip < 0) {
      // B operation
      addCigarState(-skip, ActionsHelper.CG_OVERLAP_IN_READ);
    }
  }

  @Override
  public boolean isPaired() {
    return true;
  }

  @Override
  public PrereadType prereadType() {
    return PrereadType.CG;
  }

  @Override
  public void processFragment(String id, int fragmentStart, byte[] data, int length) throws IOException {
    reseedErrorRandom(mFrameRandom.nextLong());
    final boolean forwardFrame = mFrameRandom.nextBoolean();

    final String nameLeft = generateRead(id, fragmentStart, data, length, forwardFrame, true);

    if (mReadBytesUsed == mReadBytes.length) {
      mReadWriter.writeLeftRead(nameLeft, mReadBytes, mQualityBytes, mReadBytes.length);
      mResidueCount += mReadBytes.length;
    } else {
      throw new FragmentTooSmallException(length, mReadBytes.length);
    }

    final String nameRight = generateRead(id, fragmentStart, data, length, forwardFrame, false);
    if (mReadBytesUsed == mReadBytes.length) {
      mReadWriter.writeRightRead(nameRight, mReadBytes, mQualityBytes, mReadBytes.length);
      mResidueCount += mReadBytes.length;
    } else {
      throw new FragmentTooSmallException(length, mReadBytes.length);
    }
  }

  protected abstract String generateRead(String id, int fragmentStart, byte[] data, int length, boolean forwardFrame, boolean leftArm);

}
