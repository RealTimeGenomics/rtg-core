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
package com.rtg.variant.bayes.snp;

import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;
import com.rtg.variant.bayes.ReferenceBasedFactory;

/**
 * Holds a set of objects that act as if they are mapped to an array
 * but which use a circular buffer to hold only a small number
 * of copies.
 * @param <D> type of object
 */
public class ReferenceBasedBuffer<D> extends IntegralAbstract {

  private static final int INITIAL_LENGTH = 150;

  protected final ReferenceBasedFactory<D> mFactory;

  private D[] mBuffer;

  protected final byte[] mTemplate;

  private int mBase;

  private int mCurrent;

  /**
   * @param factory for creating objects in buffer.
   * @param template nucleotides.
   * @param start of region being processed on template (0 based)
   */
  public ReferenceBasedBuffer(final ReferenceBasedFactory<D> factory, final byte[] template, int start) {
    this(INITIAL_LENGTH, factory, template, start);
  }

  /**
   * @param initialCapacity initial buffer length
   * @param factory for creating objects in buffer.
   * @param template nucleotides.
   * @param start of region being processed on template (0 based)
   */
  public ReferenceBasedBuffer(int initialCapacity, final ReferenceBasedFactory<D> factory, final byte[] template, int start) {
    mFactory = factory;
    mBuffer = makeArray(initialCapacity);
    mTemplate = template;
    mBase = start;
    mCurrent = 0;
  }

  protected D[] makeArray(int length) {
    @SuppressWarnings("unchecked")
    final D[] theArray = (D[]) new Object[length];
    return theArray;
  }

  /**
   * Get the current base position.
   * @return the current base position.
   */
  public int base() {
    return mBase;
  }

  /**
   * Get the object specified by index.
   * @param index in range of external indices, used to retrieve from buffer.
   * @return object held in buffer.
   * @throws IllegalArgumentException if index less than the current base.
   */
  public D get(final int index) {
    final int i = find(index);
    final D model = mBuffer[i];
    final D res;
    if (model != null) {
      res = model;
    } else {
      res = make(index);
      mBuffer[i] = res;
    }
    return res;
  }

  // Make a D that is appropriate to the current position on the reference
  protected D make(int index) {
    final int nt = mTemplate[index] - 1;
    assert -1 <= nt && nt < 4;
    return mFactory.make(nt);
  }

  /**
   * Step the base up by one location and return the object at the
   * (old) base position.<br> The buffer at this point is nulled.
   * @return current first object in buffer.
   */
  public D step() {
    final int i = find(mBase);
    assert i == mCurrent;
    final D theObj = mBuffer[i];
    final D res;
    if (theObj != null) {
      res = theObj;
      mBuffer[i] = null;
    } else {
      res = make(mBase);
    }
    mBase++;
    mCurrent++;
    if (mCurrent == mBuffer.length) {
      mCurrent = 0;
    }
    return res;
  }

  /**
   * Locate index in the internal index range (that is it is a valid reference into
   * our buffer). May as a side-effect resize the buffer and alter the current position.
   * @param index external index value.
   * @return the internal index.
   */
  int find(final int index) {
    final int j = index - mBase;
    if (j < 0) {
      throw new IllegalArgumentException("Index less than base. index=" + index + " base=" + mBase);
    }
    resize(j + 1);
    assert j < mBuffer.length;
    final int k = mCurrent + j;
    final int l;
    if (k >= mBuffer.length) {
      l = k - mBuffer.length;
      assert l >= 0 && l < mCurrent;
    } else {
      l = k;
    }
    assert l >= 0 && l < mBuffer.length;
    return l;
  }

  /**
   * Ensure the buffer has sufficient locations to hold at least
   * target number of objects.
   * @param targetLength minimal length of buffer after call.
   */
  private void resize(final int targetLength) {
    if (targetLength <= mBuffer.length) {
      return;
    }
    final int t = (targetLength + 1) * 2;
    final D[] newBuffer = makeArray(t);
    final int endLength = mBuffer.length - mCurrent;
    System.arraycopy(mBuffer, mCurrent, newBuffer, 0, endLength);
    System.arraycopy(mBuffer, 0, newBuffer, endLength, mCurrent);
    //    for (int i = mCurrent, j = 0; i < mBuffer.length; i++, j++) {
    //      newBuffer[j] = mBuffer[i];
    //    }
    //    for (int i = 0, j = mBuffer.length - mCurrent; i < mBuffer.length; i++, j++) {
    //      newBuffer[j] = mBuffer[i];
    //    }
    mCurrent = 0;
    mBuffer = newBuffer;
    assert globalIntegrity();
    assert targetLength <= mBuffer.length;
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("Buffer length=").append(mBuffer.length).append(" base=").append(mBase).append(" current=").append(mCurrent).append(LS);
    for (int i = 0; i < mBuffer.length; i++) {
      if (i == mCurrent) {
        sb.append(">>>>").append(LS);
      }
      sb.append("[").append(i).append("] ");
      final D model = mBuffer[i];
      sb.append("    ").append(model == null ? "null" : model.toString()).append(LS);
    }
    sb.append(LS);
  }

  @Override
  public boolean integrity() {
    Exam.assertNotNull(mFactory);
    Exam.assertNotNull(mBuffer);
    Exam.assertTrue(mBase >= 0);
    Exam.assertTrue(mCurrent >= 0 && mCurrent < mBuffer.length);
    return true;
  }
}
