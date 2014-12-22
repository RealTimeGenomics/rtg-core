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
package com.rtg.util.bytecompression;

import java.io.IOException;
import java.io.InputStream;

/**
 */
final class SingleByteArray extends ByteArray {
  private final byte[] mData;

  public SingleByteArray(final int size) {
    mData = new byte[size];
  }

  @Override
  public byte get(long offset) {
    return mData[(int) offset];
  }

  @Override
  public void get(byte[] dest, long offset, int count) {
    System.arraycopy(mData, (int) offset, dest, 0, count);
  }

  @Override
  public void set(long offset, byte value) {
    mData[(int) offset] = value;
  }

  @Override
  public void set(long offset, byte[] buffer, int count) {
    set(offset, buffer, 0, count);
  }

  @Override
  public void set(long offset, byte[] buffer, int bOffset, int count) {
    System.arraycopy(buffer, bOffset, mData, (int) offset, count);
  }

  public void load(final InputStream stream, final long offset, final int count) throws IOException {
    final int length = stream.read(mData, (int) offset, count);
    if (length != count) {
      throw new IOException();
    }
  }

  @Override
  public long length() {
    return mData.length;
  }
}
