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

import com.rtg.reader.SdfId;
import com.rtg.util.PortableRandom;

/**
 */
public class UnknownBaseReadWriter implements ReadWriter {
  final ReadWriter mInternal;
  final double mRate;
  byte[] mData;
  final PortableRandom mRandom;
  /**
   * put N in reads at a specified frequency
   * @param internal where to write the N filled reads
   * @param rate the rate to insert N
   * @param random a source of randomness
   */
  public UnknownBaseReadWriter(ReadWriter internal, double rate, PortableRandom random) {
    mInternal = internal;
    mRate = rate;
    mRandom = random;
  }

  @Override
  public void close() throws IOException {
    mInternal.close();

  }

  @Override
  public void identifyTemplateSet(SdfId... templateIds) {
    mInternal.identifyTemplateSet(templateIds);
  }

  @Override
  public void identifyOriginalReference(SdfId referenceId) {
    mInternal.identifyOriginalReference(referenceId);
  }

  private void filter(byte[] data, int length) {
    if (mData == null || length > mData.length) {
      mData = new byte[length];
    }
    for (int i = 0; i < length; i++) {
      if (mRandom.nextDouble() < mRate) {
        mData[i] = 0;
      } else {
        mData[i] = data[i];
      }
    }
  }

  @Override
  public void writeRead(String name, byte[] data, byte[] qual, int length) throws IOException {
    filter(data, length);
    mInternal.writeRead(name, mData, qual, length);
  }

  @Override
  public void writeLeftRead(String name, byte[] data, byte[] qual, int length) throws IOException {
    filter(data, length);
    mInternal.writeLeftRead(name, mData, qual, length);
  }

  @Override
  public void writeRightRead(String name, byte[] data, byte[] qual, int length) throws IOException {
    filter(data, length);
    mInternal.writeRightRead(name, mData, qual, length);
  }

  @Override
  public int readsWritten() {
    return mInternal.readsWritten();
  }

}
