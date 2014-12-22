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

package com.rtg.util.iterators;

import java.util.Iterator;


/**
 * Given an <code>Iterator&lt;X&gt;</code> and a transformation method <code>subIterator</code> that constructs an
 * <code>Iterator&lt;Y&gt;</code> from an instance of <code>X</code> subclasses will act as an iterator of all
 * <code>Y</code> that can be produced
 * @param <X> type of source objects
 * @param <Y> resulting type objects
 */
public final class ComposeIterators<X, Y> implements Iterator<Y> {

  private final Iterator<X> mSourceIterator;
  private final Transform<X, Iterator<Y>> mTransform;
  private Iterator<Y> mCurrent;

  ComposeIterators(Iterator<X> sourceIterator, Transform<X, Iterator<Y>> transform) {
    mSourceIterator = sourceIterator;
    mTransform = transform;
  }

  private void stepWhile() {
    while (mCurrent != null && !mCurrent.hasNext()) {
      step();
    }
  }

  private void step() {
    if (mSourceIterator.hasNext()) {
      mCurrent = mTransform.trans(mSourceIterator.next());
    } else {
      mCurrent = null;
    }
  }

  @Override
  public boolean hasNext() {
    if (mCurrent == null) {
      step();
      stepWhile();
    }
    if (mCurrent == null) {
      return false;
    } else {
      return mCurrent.hasNext();
    }
  }

  @Override
  public Y next() {
    final Y res = mCurrent.next();
    stepWhile();
    return res;
  }

  @Override
  public void remove() {
    throw new UnsupportedOperationException();
  }
}
