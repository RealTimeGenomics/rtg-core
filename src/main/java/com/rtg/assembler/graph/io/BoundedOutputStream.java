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

package com.rtg.assembler.graph.io;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.OutputStream;

/**
 */
public abstract class BoundedOutputStream extends OutputStream {

  private final ByteArrayOutputStream mBuffer = new ByteArrayOutputStream();

  private final int mMaxLength;

  private OutputStream mProxy = null;

  private int mTotalSize = 0;

  private int mCount = 1;

  BoundedOutputStream(final int maxLength) {
    mMaxLength = maxLength;
  }

  abstract OutputStream nextProxy(final int count) throws IOException;

  @Override
  public void close() throws IOException {
    flush();
    mProxy.close();
  }

  @Override
  public void flush() throws IOException {
    final int size = mBuffer.size();
    if (size > mMaxLength) {
      throw new RuntimeException(size + ":" + mMaxLength);
    }
    final long nextSize = mTotalSize + size;
    //System.err.println("nextSize=" + nextSize + " mTotalSize=" + mTotalSize + " size=" + size + " mMaxLength=" + mMaxLength);
    if (nextSize > mMaxLength) {
      mProxy.close();
      mProxy = null;
      mTotalSize = 0;
    }
    if (mProxy == null) {
      mProxy = nextProxy(mCount);
      ++mCount;
    }
    mProxy.write(mBuffer.toByteArray());
    mTotalSize += size;
    mProxy.flush();
    mBuffer.reset();
  }

  @Override
  public void write(byte[] b, int off, int len) {
    mBuffer.write(b, off, len);
  }


  @Override
  public void write(byte[] b) throws IOException {
    mBuffer.write(b);
  }

  @Override
  public void write(int arg0) {
    mBuffer.write(arg0);
  }

}
