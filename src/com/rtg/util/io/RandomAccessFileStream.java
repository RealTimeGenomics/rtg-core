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

package com.rtg.util.io;

import java.io.IOException;
import java.io.RandomAccessFile;

/**
 * SeekableStream based on random access file
 */
public class RandomAccessFileStream extends SeekableStream {

  private final RandomAccessFile mRaf;

  /**
   * Constructor
   * @param raf base random access file
   */
  public RandomAccessFileStream(RandomAccessFile raf) {
    mRaf = raf;
  }

  @Override
  public long getPosition() throws IOException {
    return mRaf.getFilePointer();
  }

  @Override
  public long length() throws IOException {
    return mRaf.length();
  }

  @Override
  public void seek(long pos) throws IOException {
    mRaf.seek(pos);
  }

  @Override
  public int read() throws IOException {
    return mRaf.read();
  }

  @Override
  public int available() throws IOException {
    final long ret = mRaf.length() - mRaf.getFilePointer();
    if (ret > Integer.MAX_VALUE) {
      return Integer.MAX_VALUE;
    }
    return (int) ret;
  }

  /**
   */
  @Override
  public void close() throws IOException {
    mRaf.close();
  }

  @Override
  public synchronized void mark(int readlimit) {
    throw new UnsupportedOperationException();
  }

  @Override
  public boolean markSupported() {
    return false;
  }

  @Override
  public int read(byte[] b) throws IOException {
    return mRaf.read(b);
  }

  @Override
  public int read(byte[] b, int off, int len) throws IOException {
    return mRaf.read(b, off, len);
  }

  @Override
  public synchronized void reset() {
    throw new UnsupportedOperationException();
  }

  @Override
  public long skip(long n) throws IOException {
    final long skip = n > available() ? (long) available() : n;
    mRaf.seek(mRaf.getFilePointer() + skip);
    return skip;
  }

}
