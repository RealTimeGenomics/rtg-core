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
package com.rtg.index.hash.ngs.protein;

import java.io.IOException;

import com.rtg.index.hash.ngs.ReadCall;
import com.rtg.index.hash.ngs.general.SingleMask;
import com.rtg.util.integrity.Exam;

/**
 */
public final class ProteinExtractRead extends ProteinExtractAbstract {

  private final ReadCall mReadCall;

  private final int mIndex;

  private int mReadId = -1;

  ProteinExtractRead(final SingleMask singleMask, final ReadCall readCall, final int index) {
    super(singleMask);
    mReadCall = readCall;
    mIndex = index;
    integrity();
  }

  void readCall(int id, long[] vs) throws IOException {
    mReadId = id;
    //this ends up calling masked() below and uses the mReadId
    mask(vs);
  }

  @Override
  protected void masked(final long bits) {
    if (bits == 0) { //if resulting hash is zero do not add into index, in protein in translate to all X's
      return;
    }
    mReadCall.readCall(mReadId, bits, mIndex);
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("ProteinExtractRead readCall=").append(mReadCall).append(" index=").append(mIndex).append(" id=").append(mReadId);
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mReadCall != null);
    Exam.assertTrue(mIndex >= 0);
    return true;
  }
}
