/*
 * Copyright (c) 2016. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */

package com.rtg.reader;

import java.util.Iterator;

import htsjdk.samtools.util.RuntimeIOException;

/**
 */
public class FastqPairIterator implements Iterator<FastqPair> {
  private final Iterator<FastqSequence> mLeft;
  private final Iterator<FastqSequence> mRight;

  FastqPairIterator(Iterator<FastqSequence> left, Iterator<FastqSequence> right) {
    mLeft = left;
    mRight = right;
  }

  @Override
  public boolean hasNext() {
    if (mLeft.hasNext() != mRight.hasNext()) {
      throw new RuntimeIOException("Left reader has different number of sequences from right reader");
    }
    return  mLeft.hasNext();
  }

  @Override
  public FastqPair next() {
    return new FastqPair(mLeft.next(), mRight.next());
  }
}
