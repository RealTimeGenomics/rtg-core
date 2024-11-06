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
