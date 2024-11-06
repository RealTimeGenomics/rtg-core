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
package com.rtg.ngs.tempstage;

import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.Channels;
import java.nio.channels.ReadableByteChannel;

import com.reeltwo.jumble.annotations.TestClass;

/**
 */
@TestClass("com.rtg.ngs.tempstage.TempRecordNioTest")
public class TempRecordReaderNio implements TempRecordReader {

  private final ReadableByteChannel mChannel;
  private final RecordFactory mFact;
  private final ByteBuffer mBuffer;

  /**
   * @param inputStream stream which will be read from
   * @param fact factory to create temp file records
   */
  public TempRecordReaderNio(InputStream inputStream, RecordFactory fact) {
    mChannel = Channels.newChannel(inputStream);
    mFact = fact;
    mBuffer = ByteBuffer.allocate(64 * 1024);
    mBuffer.order(ByteOrder.nativeOrder());
    mBuffer.flip();
  }

  @Override
  public BinaryTempFileRecord readRecord() throws IOException {
    final BinaryTempFileRecord ret = mFact.createRecord();
    ret.readNio(mBuffer, mChannel);
    if (!ret.isSentinelRecord()) {
      return ret;
    }
    return null;
  }

  @Override
  public void close() throws IOException {
    mChannel.close();
  }
}
