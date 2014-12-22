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

import com.rtg.index.hash.ngs.TemplateCall;
import com.rtg.index.hash.ngs.general.SingleMask;
import com.rtg.util.integrity.Exam;

/**
 */
public final class ProteinExtractTemplate extends ProteinExtractAbstract {

  final TemplateCall mTemplateCall;

  private final int mIndex;

  private int mEndPosition = -1;

  ProteinExtractTemplate(final SingleMask singleMask, final TemplateCall templateCall, final int index) {
    super(singleMask);
    mTemplateCall = templateCall;
    mIndex = index;
    integrity();
  }

  void templateCall(int endPosition, long[] vs) throws IOException {
    mEndPosition = endPosition;
    //this ends up calling masked() below and uses the mReadId
    maskIndel(vs);
  }

  @Override
  protected void masked(final long bits) throws IOException {
    if (bits == 0) { //if resulting hash is zero do not add into index, in protein in translate to all X's
      return;
    }
    //System.err.println("hash=" + com.rtg.util.Utils.toBitsSep(bits) + " " + bits);
    mTemplateCall.templateCall(mEndPosition, bits, mIndex);
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("ProteinExtractTemplate templateCall=").append(mTemplateCall).append(" index=").append(mIndex).append(" end=").append(mEndPosition);
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mTemplateCall != null);
    Exam.assertTrue(mIndex >= 0);
    return true;
  }
}
