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

package com.rtg.assembler;

import java.util.Iterator;

import com.rtg.mode.DNA;

/**
*/
final class KmerIterator implements Iterator<Kmer> {

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
    ++mStart;
    while (mI - mStart < mKmerSize && hasNext()) {
      ++mI;
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
