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
import java.util.Collection;

import com.rtg.index.hash.ngs.HashFunctionFactory;
import com.rtg.index.hash.ngs.ImplementHashFunction;
import com.rtg.index.hash.ngs.NgsHashFunction;
import com.rtg.index.hash.ngs.ReadCall;
import com.rtg.index.hash.ngs.TemplateCall;
import com.rtg.index.hash.ngs.general.SingleMask;
import com.rtg.index.hash.ngs.general.Skeleton;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.integrity.Exam;

/**
 * General mask that does protein matching using a given skeleton.
 * A set of all necessary skeletons are computed from the -a, -b, -c parameters.
 *
 */
public class ProteinMask extends ImplementHashFunction {

  /** Number of bits to store one amino acid */
  static final int PROTEIN_BITS = 5;

  private static final class ProteinMaskFactory implements HashFunctionFactory {

    private final Skeleton mSkeleton;

    ProteinMaskFactory(Skeleton skeleton) {
      mSkeleton = skeleton;
    }

    @Override
    public NgsHashFunction create(final ReadCall readCall, final TemplateCall templateCall) {
      //System.err.println("create " + mSkeleton);
      final ProteinMask mask = new ProteinMask(mSkeleton, readCall, templateCall);
      assert mask.integrity();
      return mask;
    }

    @Override
    public int hashBits() {
      return mSkeleton.windowLength() * PROTEIN_BITS;
    }

    @Override
    public int numberWindows() {
      return mSkeleton.numberWindows();
    }

    @Override
    public int windowBits() {
      return mSkeleton.windowLength() * PROTEIN_BITS;
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
    return new ProteinMaskFactory(sk);
  }

  private final Skeleton mSkeleton;

  private final int mNumberMasks;

  private ProteinExtractRead[] mReadMasks;

  private ProteinExtractTemplate[] mTemplateMasks;

  /**
   * Encodes up to 64 amino acids, from most significant bits to least significant.
   * Eg. the amino acid Tryptophan (W) has ordinal value 19 = 16+0+0+2+1,
   * so if this is stored in bit B (some power of 2), then
   * <pre>
   *  mValuesF[0]&B != 0 // for the 16
   *  mValuesF[1]&B == 0 // for the 8
   *  mValuesF[2]&B == 0 // for the 4
   *  mValuesF[3]&B != 0 // for the 2
   *  mValuesF[4]&B != 0 // for the 1
   * </pre>
   */
  private long[] mValuesF = new long[PROTEIN_BITS];

  @Override
  public ProteinMask clone() throws CloneNotSupportedException {
    final ProteinMask clone = (ProteinMask) super.clone();
    clone.mValuesF = mValuesF.clone();
    return clone;
  }

  /**
   * @param sk the skeleton used to generate the mask.
   * @param readCall used in subclasses to process results of read hits.
   * @param templateCall used in subclasses to process results of template hits.
   */
  protected ProteinMask(final Skeleton sk, final ReadCall readCall, final TemplateCall templateCall) {
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
    mReadMasks = new ProteinExtractRead[mNumberMasks];
    for (int i = 0; i < mNumberMasks; ++i) {
      mReadMasks[i] = new ProteinExtractRead(ma[i], mReadCall, i);
    }
    //System.err.println("set masks HashFunction:" + System.identityHashCode(this) + " mTemplateCall:" + System.identityHashCode(mTemplateCall));
    mTemplateMasks = new ProteinExtractTemplate[mNumberMasks];
    for (int i = 0; i < mNumberMasks; ++i) {
      mTemplateMasks[i] = new ProteinExtractTemplate(ma[i], mTemplateCall, i);
    }
  }

  @Override
  public void templateForward(int endPosition) throws IOException {
    if (!isValid()) {
      //System.err.println("WE return");
      return;
    }
    mTemplateCall.setReverse(false);
    templateAll(endPosition, mValuesF);
  }

  /**
   * Protein does not support reverse frame.
   */
  @Override
  public void templateReverse(int endPosition) {
    throw new UnsupportedOperationException();
  }

  @Override
  public void readAll(int readId, boolean reverse) throws IOException {
    //System.err.println("readAll readId=" + readId + " reverse=" + reverse);
    assert !reverse;
    readAll(readId, mValuesF);
  }


  /**
   * Similar to the usual <code>templateAll</code>, but for proteins.
   * @param endPosition end position
   * @param vs an array of 5 longs.
   * @throws IOException If an I/O error occurs
   */
  public void templateAll(int endPosition, long[] vs) throws IOException {
    //System.err.println("templateAll endPosition=" + endPosition);
    if (!isValid()) { //mSoFar < readLength()) {
      //System.err.println("isValid=" + true);
      return;
    }
    for (int i = 0; i < mNumberMasks; ++i) {
      mTemplateMasks[i].templateCall(endPosition, vs);
    }
    mTemplateCall.done();
  }

  /**
   * Return whether this mask can successfully <code>readAll</code> or <code>templateAll</code>
   * @return true for yes
   */
  public boolean isValid() {
    return mSoFar >= mSkeleton.windowLength();
    //    return mSoFar >= readLength();
  }


  /**
   * Similar to the usual <code>readAll</code>, but for proteins.
   * @param readId read id
   * @param vs an array of 5 longs.
   * @throws IOException If an I/O error occurs
   */
  public void readAll(final int readId, final long[] vs) throws IOException {
    if (!isValid()) { // mSoFar < readLength()) {
      return;
    }
    for (int i = 0; i < mNumberMasks; ++i) {
      mReadMasks[i].readCall(readId, vs);
    }
  }

  @Override
  public void hashStep(final byte code) {
    int fwd = code;
    assert 0 <= code && code < (1 << PROTEIN_BITS);
    for (int i = PROTEIN_BITS - 1; i >= 0; --i) {
      mValuesF[i] = (mValuesF[i] << 1) | (fwd & 1);
      fwd = fwd >>> 1;
    }

    ++mSoFar;
//    System.err.println("hashStep code=" + code);
//    for (int i = 0; i < PROTEIN_BITS; ++i) {
//      System.err.println("mValuesF[" + i + "]=" + com.rtg.util.Utils.toBitsSep(mValuesF[i]));
//    }
  }

  @Override
  public int numberWindows() {
    return mNumberMasks;
  }

  @Override
  public void toString(StringBuilder sb) {
    sb.append("Mask l=").append(readLength()).append(" w=").append(mSkeleton.windowBits() / 2).append(" s=").append(mSkeleton.substitutions()).append(" i=").append(mSkeleton.indels());
  }

  boolean integrity(int soFar, long... vs) {
    Exam.assertEquals(soFar, mSoFar);
    Exam.assertEquals(vs.length, mValuesF.length);
    for (int i = 0; i < mValuesF.length; ++i) {
      Exam.assertEquals(vs[i], mValuesF[i]);
    }
    return true;
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

  @Override
  public void reset() {
    mSoFar = 0;
    for (int i = 0; i < PROTEIN_BITS; ++i) {
      mValuesF[i] = 0L;
    }
  }

  /**
   * Not supported for proteins.
   */
  @Override
  public void setValues(final int id2, boolean reverse) {
    throw new UnsupportedOperationException();
  }

  /**
   * Not supported for proteins.
   */
  @Override
  public int fastScore(final int readId) {
    return 0;
    //throw new UnsupportedOperationException();
  }

  /**
   * Not supported for proteins.
   */
  @Override
  public int indelScore(final int readId) {
    return 0;
    //throw new UnsupportedOperationException();
  }

  /**
   * Not supported for proteins.
   */
  @Override
  public void readAll(int readId, long v0, long v1) {
    throw new UnsupportedOperationException();
  }

  /**
   * Not supported for proteins.
   */
  @Override
  public void templateAll(int endPosition, long v0, long v1) {
    throw new UnsupportedOperationException();
  }
}
