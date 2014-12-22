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

import com.reeltwo.jumble.annotations.TestClass;


/**
 * @param <X> the class being transformed from.
 * @param <Y> the class being transformed to.
 */
@TestClass({"com.rtg.util.iterators.TransformTest", "com.rtg.util.iterators.ArrayToIteratorTest", "com.rtg.util.iterators.ComposeIteratorsTest"})
public abstract class Transform<X, Y> {

  /**
   * @param x value to be transformed.
   * @return value after transformation.
   */
  public abstract Y trans(X x);

  private static final class IterTrans<X, Y> implements Iterator<Y> {
    private final Transform<X, Y> mTrans;
    private final Iterator<X> mIt;

    /**
     * @param trans transform to be applied to the iterator.
     * @param it iterator to be transformed.
     */
    IterTrans(final Transform<X, Y> trans, final Iterator<X> it) {
      mTrans = trans;
      mIt = it;
    }

    @Override
    public boolean hasNext() {
      return mIt.hasNext();
    }

    @Override
    public Y next() {
      return mTrans.trans(mIt.next());
    }

    @Override
    public void remove() {
      mIt.remove();
    }
  }

  /**
   * Transform an iterator.
   * @param xit iterator to be transformed.
   * @return an iterator all the objects of which have been transformed.
   */
  public Iterator<Y> trans(Iterator<X> xit) {
    return new IterTrans<>(this, xit);
  }

  private static final class TransCompose<X, Y, Z> extends Transform<X, Z> {

    private final Transform<X, Y> mXY;
    private final Transform<Y, Z> mYZ;
    TransCompose(Transform<X, Y> xy, Transform<Y, Z> yz) {
      mXY = xy;
      mYZ = yz;
    }
    @Override
    public Z trans(X x) {
      return mYZ.trans(mXY.trans(x));
    }
  }

  /**
   * Compose two transforms to form a single new transform.
   * @param xy transform from x to y.
   * @param yz transform from y to z.
   * @return a transform from x to z.
   * @param <X> type of original values.
   * @param <Y> type of intermediate values.
   * @param <Z> type of final values.
   */
  public static <X, Y, Z> Transform<X, Z> compose(Transform<X, Y> xy, Transform<Y, Z> yz) {
    return new TransCompose<>(xy, yz);
  }

  /**
   * Given an iterator and a transform that takes each value in the iterator to another iterator,
   * flatten all this to give a single iterator over the final values.
   * @param x the initial iterator.
   * @param xiy the transformer that takes initial values to iterators over the final values.
   * @return a single iterator over the final values.
   * @param <X> the type of the initial values.
   * @param <Y> the type of the final values.
   */
  public static <X, Y> Iterator<Y> flatten(Iterator<X> x, Transform<X, Iterator<Y>> xiy) {
    return new ComposeIterators<>(x, xiy);
  }

  /**
   * Given an iterator over iterables over values,
   * flatten all this to give a single iterator over the values.
   * @param x the iterator over <code>Iterables</code>.
   * @return a single iterator over the final values.
   * @param <I> the iterable type.
   * @param <Y> the type of the final values.
   */
  public static <I extends Iterable<Y>, Y> Iterator<Y> flatten(Iterator<I> x) {
    final Transform<I, Iterator<Y>> trans = new Transform<I, Iterator<Y>>() {
      @Override
      public Iterator<Y> trans(I x) {
        return x.iterator();
      }
    };
    return new ComposeIterators<>(x, trans);
  }

  /**
   * Convert an array to an iterator over those values.
   * This is handy for testing.
   * @param array the array of values.
   * @return an iterator over the values.
   * @param <X> the type of the values.
   */
  public static <X> Iterator<X> array2Iterator(X[] array) {
    return new ArrayToIterator<>(array);
  }
}
