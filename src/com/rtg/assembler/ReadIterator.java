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
import java.util.List;
import java.util.NoSuchElementException;

/**
 *  Creates an iterator from a sequences reader.
 *  Will throw runtime exceptions in the event of IO failure. Ugly but I blame the Iterator interface.
 */
class ReadIterator implements Iterator<List<byte[]>> {
  private final AsyncReadSource mReader;
  private List<byte[]> mFragments;

  ReadIterator(AsyncReadSource reader) {
    mReader = reader;
    mFragments =  mReader.nextFragments();
  }


  @Override
  public boolean hasNext() {
    return mFragments != null;
  }

  @Override
  public List<byte[]> next() {
    if (mFragments == null) {
      throw new NoSuchElementException();
    }
    final List<byte[]> fragments = mFragments;
    mFragments = mReader.nextFragments();
    return fragments;
  }

  @Override
  public void remove() {
    throw new UnsupportedOperationException();
  }
}
