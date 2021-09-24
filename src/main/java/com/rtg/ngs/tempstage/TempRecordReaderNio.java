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
