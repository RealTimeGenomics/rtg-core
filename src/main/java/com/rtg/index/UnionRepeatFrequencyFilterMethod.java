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
