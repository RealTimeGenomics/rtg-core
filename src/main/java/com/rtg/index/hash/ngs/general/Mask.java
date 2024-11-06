/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
package com.rtg.index.hash.ngs.general;

import java.io.IOException;
import java.util.Collection;

import com.rtg.index.hash.ngs.HashFunctionFactory;
import com.rtg.index.hash.ngs.ImplementHashFunction;
import com.rtg.index.hash.ngs.NgsHashFunction;
import com.rtg.index.hash.ngs.ReadCall;
import com.rtg.index.hash.ngs.TemplateCall;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.integrity.Exam;

/**
 * General mask that does matching using a given skeleton.
 * A set of all necessary skeletons are computed from the -a, -b, -c parameters.
 *
 */
public final class Mask extends ImplementHashFunction {

  private static final class MaskFactory implements HashFunctionFactory {

    private final Skeleton mSkeleton;

    MaskFactory(final Skeleton skeleton) {
      mSkeleton = skeleton;
    }

    @Override
    public NgsHashFunction create(final ReadCall readCall, final TemplateCall templateCall) {
      //System.err.println("create " + mSkeleton);
      final Mask mask = new Mask(mSkeleton, readCall, templateCall);
      assert mask.integrity();
      return mask;
    }

    @Override
    public int hashBits() {
      return mSkeleton.windowBits();
    }

    @Override
    public int numberWindows() {
      return mSkeleton.numberWindows();
    }

    @Override
    public int windowBits() {
      return mSkeleton.windowBits();
    }

    /**
     * Return the window size
     * @return the size
     */
    @Override
    public int windowSize() {
      return mSkeleton.windowLength();
    }

  }

  /**
   * Construct a factory based on the skeleton.
   * @param sk to be used for factory and ultimately for the mask.
   * @return factory.
   */
  public static HashFunctionFactory factory(final Skeleton sk) {
    if (!sk.valid()) {
      throw new RuntimeException("Invalid skeleton:" + StringUtils.LS + sk);
    }
    Diagnostic.developerLog(sk.dumpMask());
    //System.err.println("skeleton=" + sk);
    return new MaskFactory(sk);
  }

  private final Skeleton mSkeleton;

  private final int mNumberMasks;

  private ExtractRead[] mReadMasks;

  private ExtractTemplate[] mTemplateMasks;

  /**
   * @param sk the skeleton used to generate the mask.
   * @param readCall used in subclasses to process results of read hits.
   * @param templateCall used in subclasses to process results of template hits.
   */
  private Mask(final Skeleton sk, final ReadCall readCall, final TemplateCall templateCall) {
    super(sk.readLength(), sk.windowLength(), readCall, templateCall);

    mSkeleton = sk;
    mNumberMasks = mSkeleton.numberWindows();
    setHashFunction();
    assert integrity();
  }


  @Override
  protected void setHashFunction() {
    super.setHashFunction();
    setMasks();
  }

  private void setMasks() {
    final Collection<SingleMask> masks = mSkeleton.masks();
    //System.err.println(masks);
    assert masks.size() == mNumberMasks;
    final SingleMask[] ma = masks.toArray(new SingleMask[mNumberMasks]);
    mReadMasks = new ExtractRead[mNumberMasks];
    for (int i = 0; i < mNumberMasks; ++i) {
      mReadMasks[i] = new ExtractRead(ma[i], mReadCall, i);
    }
    //System.err.println("set masks HashFunction:" + System.identityHashCode(this) + " mTemplateCall:" + System.identityHashCode(mTemplateCall));
    mTemplateMasks = new ExtractTemplate[mNumberMasks];
    for (int i = 0; i < mNumberMasks; ++i) {
      mTemplateMasks[i] = new ExtractTemplate(ma[i], mTemplateCall, i);
    }
  }

  @Override
  public void readAll(final int readId, long v0, long v1) throws IOException {
    if (mSoFar < readLength()) {
      return;
    }

    for (int i = 0; i < mNumberMasks; ++i) {
      //System.err.println("readAll mask i=" + i);
      mReadMasks[i].readCall(readId, v0, v1);
    }
  }

  @Override
  public void templateAll(final int endPosition, long v0, long v1) throws IOException {
    //System.err.println("templateAll endPosition=" + endPosition);
    if (mSoFar < readLength()) {
      return;
    }

    for (int i = 0; i < mNumberMasks; ++i) {
      //System.err.println("templateAll mask i=" + i);
      mTemplateMasks[i].templateCall(endPosition, v0, v1);
    }

    mTemplateCall.done();
}

  @Override
  public int numberWindows() {
    return mNumberMasks;
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("Mask l=").append(readLength()).append(" w=").append(mSkeleton.windowBits() / 2).append(" s=").append(mSkeleton.substitutions()).append(" i=").append(mSkeleton.indels());
  }



  @Override
  public boolean integrity() {
    super.integrity();
    Exam.assertTrue(mNumberMasks > 0);
    Exam.assertTrue(mReadMasks != null && mReadMasks.length == mNumberMasks);
    Exam.assertTrue(mTemplateMasks == null ? "null" : "length=" + mTemplateMasks.length + " numberMasks=" + mNumberMasks, mTemplateMasks != null && mTemplateMasks.length == mNumberMasks);
    for (int i = 0; i < mNumberMasks; ++i) {
      Exam.assertTrue(mReadMasks[i] != null);
      Exam.assertTrue(mTemplateMasks[i] != null && mTemplateMasks[i].mTemplateCall == mTemplateCall);
    }
    return true;
  }
}
