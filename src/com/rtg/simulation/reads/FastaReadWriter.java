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

import com.rtg.mode.DnaUtils;
import com.rtg.reader.SdfId;

/**
 * FASTA style read simulator output
 */
public class FastaReadWriter implements ReadWriter {

  private final Appendable mAppend;
  private int mTotal = 0;
  private boolean mExpectLeft = true;

  /**
   * Constructor
   * @param append destination for output
   */
  public FastaReadWriter(Appendable append) {
    mAppend = append;
  }

  @Override
  public void identifyTemplateSet(SdfId... templateIds) {
    // Ignored
  }

  @Override
  public void identifyOriginalReference(SdfId referenceId) {
    // Ignored
  }

  @Override
  public void writeRead(String name, byte[] data, byte[] qual, int length) throws IOException {
    writeSequence(mTotal + " " + name, data, length);
    ++mTotal;
  }

  @Override
  public void writeLeftRead(String name, byte[] data, byte[] qual, int length) throws IOException {
    if (!mExpectLeft) {
      throw new IllegalStateException();
    }
    writeSequence(mTotal + " " + name + "/Left", data, length);
    mExpectLeft = !mExpectLeft;
  }

  @Override
  public void writeRightRead(String name, byte[] data, byte[] qual, int length) throws IOException {
    if (mExpectLeft) {
      throw new IllegalStateException();
    }
    writeSequence(mTotal + " " + name + "/Right", data, length);
    mExpectLeft = !mExpectLeft;
    ++mTotal;
  }

  private void writeSequence(String name, byte[] data, int length) throws IOException {
    mAppend.append(">");
    mAppend.append(name);
    mAppend.append("\n");
    mAppend.append(DnaUtils.bytesToSequenceIncCG(data, 0, length));
    mAppend.append("\n");
  }

  @Override
  public void close() {
  }


  @Override
  public int readsWritten() {
    return mTotal;
  }

}
