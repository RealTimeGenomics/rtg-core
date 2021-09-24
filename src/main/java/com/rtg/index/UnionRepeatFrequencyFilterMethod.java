/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */
package com.rtg.index;

import java.util.Arrays;
import java.util.Collection;

/**
 * Keep a hash only if all delegate methods say it is OK to keep.
 */
public class UnionRepeatFrequencyFilterMethod implements IndexFilterMethod {

  /**
   * Hashes with more than this number of hits are assumed to be repeats and are
   * not reported during searches. This can give performance improvements by
   * eliminating the construction of index hit objects that only contain these
   * high frequency chunks.
   */
  protected IndexFilterMethod[] mMethods;

  /**
   * @param methods the delegate filter methods
   */
  public UnionRepeatFrequencyFilterMethod(Collection<IndexFilterMethod> methods) {
    this(methods.toArray(new IndexFilterMethod[0]));
  }

  /**
   * @param methods the delegate filter methods
   */
  public UnionRepeatFrequencyFilterMethod(IndexFilterMethod... methods) {
    mMethods = methods;
  }

  @Override
  public IndexFilterMethod threadClone() {
    final IndexFilterMethod[] cloned = new IndexFilterMethod[mMethods.length];
    for (int i = 0; i < cloned.length; i++) {
      cloned[i] = mMethods[i].threadClone();
    }
    return new UnionRepeatFrequencyFilterMethod(cloned);
  }

  @Override
  public void initialize(Index index) {
    for (IndexFilterMethod m : mMethods) {
      m.initialize(index);
    }
  }

  @Override
  public boolean keepHash(long hash, long numHits) {
    for (IndexFilterMethod m : mMethods) {
      if (!m.keepHash(hash, numHits)) {
        return false;
      }
    }
    return true;
  }

  @Override
  public String toString() {
    return "Union " + Arrays.toString(mMethods);
  }
}
