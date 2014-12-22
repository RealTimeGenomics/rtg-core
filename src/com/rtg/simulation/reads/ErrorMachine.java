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
import com.rtg.reader.SdfId;
import com.rtg.util.PortableRandom;

/**
 * A read simulation machine which introduces additional errors
 * Specifically PCR amplification and Chimeras. Will not duplicate AND chimera a fragment at the same time.
 */
public class ErrorMachine implements Machine {

  private final Machine mParentMachine;
  private final PortableRandom mErrorRandom;

  private final double mPcrDupRate;
  private final double mChimeraRate;

  private int mChimeraCount = 0;
  private String mChimeraId;
  private byte[] mChimeraData;
  private byte[] mChimeraDataBuf;
  private int mChimeraLength;

  /**
   * @param randomSeed seed for random number generator
   * @param parent the parent machine to wrap
   * @param pcrDuplicationRate the rate of PCR duplication errors
   * @param chimeraRate the rate of chimeric fragment errors
   */
  public ErrorMachine(long randomSeed, Machine parent, double pcrDuplicationRate, double chimeraRate) {
    mErrorRandom = new PortableRandom(randomSeed);
    mParentMachine = parent;

    mPcrDupRate = pcrDuplicationRate;
    mChimeraRate = pcrDuplicationRate + chimeraRate;
  }

  @Override
  public void setQualRange(byte minq, byte maxq) {
    mParentMachine.setQualRange(minq, maxq);
  }

  @Override
  public void setReadWriter(ReadWriter rw) {
    mParentMachine.setReadWriter(rw);
  }

  @Override
  public void identifyTemplateSet(SdfId... templateIds) {
    mParentMachine.identifyTemplateSet(templateIds);
  }

  @Override
  public void identifyOriginalReference(SdfId referenceId) {
    mParentMachine.identifyOriginalReference(referenceId);
  }

  @Override
  public void processFragment(String id, int fragmentStart, byte[] data, int length) throws IOException {

    if (mChimeraId != null) {
      final int catLength = mChimeraLength + length;
      if (mChimeraDataBuf == null || mChimeraDataBuf.length < catLength) {
        mChimeraDataBuf = new byte[catLength];
      }
      System.arraycopy(mChimeraData, 0, mChimeraDataBuf, 0, mChimeraLength);
      System.arraycopy(data, 0, mChimeraDataBuf, mChimeraLength, length);

      mParentMachine.processFragment("chimera" + mChimeraCount + '/', Integer.MIN_VALUE, mChimeraDataBuf, catLength);

      mChimeraCount++;
      mChimeraId = null;
      return;
    }

    final double nextRandom = mErrorRandom.nextDouble();
    if (nextRandom < mPcrDupRate) {
      mParentMachine.processFragment(id, fragmentStart, data, length);
      mParentMachine.processFragment("dupe-" + id, fragmentStart, data, length); //TODO chance to get higher levels of duplication
    } else if (nextRandom < mChimeraRate) {
      //store a fragment to mangle together with the next fragment
      mChimeraId = id;
      mChimeraData = data;
      mChimeraLength = length;
    } else {
      mParentMachine.processFragment(id, fragmentStart, data, length);
    }
  }

  @Override
  public long residues() {
    return mParentMachine.residues();
  }

  @Override
  public boolean isPaired() {
    return mParentMachine.isPaired();
  }

  @Override
  public PrereadType machineType() {
    return mParentMachine.machineType();
  }

  @Override
  public String formatActionsHistogram() {
    return mParentMachine.formatActionsHistogram();
  }
}
