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
package com.rtg.index.hash.ngs;

import java.io.IOException;

import com.rtg.index.hash.ngs.general.Skeleton;
import com.rtg.index.hash.ngs.protein.ProteinMask;
import com.rtg.launcher.HashingRegion;

class FakeProteinMask extends ProteinMask {

  int mReadCalls;
  int mTemplateCalls;
  int mClones;

  public FakeProteinMask(Skeleton sk, ReadCall readCall, TemplateCall templateCall) {
    super(sk, readCall, templateCall);
  }

  @Override
  public void readAll(int readId, boolean reverse) {
    junit.framework.Assert.assertFalse(reverse);
    mReadCalls++;
  }

  @Override
  public void templateForward(int endPosition) {
    mTemplateCalls++;
  }

  @Override
  public NgsHashFunction threadClone(HashingRegion region) throws IOException {
    mClones++;
    return super.threadClone(region);
  }
}
