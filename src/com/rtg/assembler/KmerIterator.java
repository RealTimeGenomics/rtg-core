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

package com.rtg.assembler;

import java.util.Iterator;

import com.rtg.mode.DNA;

/**
*/
class KmerIterator implements Iterator<Kmer> {

  private final byte[] mRead;
  private int mI = 0;
  private int mStart = -1;
  private final int mKmerSize;
  private final KmerFactory mKmerFactory;

  KmerIterator(byte[] read, KmerFactory factory, int kmerSize) {
    mRead = read;
    mKmerFactory = factory;
    mKmerSize = kmerSize;
    step();
  }

  @Override
  public boolean hasNext() {
    return mI <= mRead.length;
  }

  private void step() {
    mStart++;
    while (hasNext() && mI - mStart < mKmerSize) {
      mI++;
      if (mI > mRead.length) {
        break;
      }
      if (mRead[mI - 1] == DNA.N.ordinal()) {
        mStart = mI;
      }
    }
  }
  @Override
  public Kmer next() {
    final Kmer kmer = mKmerFactory.make(mRead, mStart, mI);
    step();
    return kmer;
  }

  @Override
  public void remove() {
    throw new UnsupportedOperationException();
  }
}
