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
import java.io.OutputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.Channels;
import java.nio.channels.WritableByteChannel;

import com.reeltwo.jumble.annotations.TestClass;

/**
 */
@TestClass("com.rtg.ngs.tempstage.TempRecordNioTest")
public class TempRecordWriterNio implements TempRecordWriter {

  private final WritableByteChannel mOutChannel;
  private final ByteBuffer mBuffer;

  /**
   * @param out the output stream which will be written to
   */
  public TempRecordWriterNio(OutputStream out) {
    mOutChannel = Channels.newChannel(out);
    mBuffer = ByteBuffer.allocate(64 * 1024);
    mBuffer.order(ByteOrder.nativeOrder());
  }

  @Override
  public void writeRecord(BinaryTempFileRecord rec) throws IOException {
    mBuffer.clear();
    rec.writeNio(mBuffer);
    mBuffer.flip();
    mOutChannel.write(mBuffer);
  }

  @Override
  @SuppressWarnings("try")
  public void close() throws IOException {
    try (WritableByteChannel ignored = mOutChannel) {
      final BinaryTempFileRecord sent = new BinaryTempFileRecord(false, false, false, false);
      sent.setSentinelRecord();
      writeRecord(sent);
    }
  }
}
