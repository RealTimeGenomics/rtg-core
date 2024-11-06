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

import java.io.File;
import java.io.IOException;

import com.rtg.util.IORunnable;
import com.rtg.util.SimpleThreadPool;
import com.rtg.util.io.FileUtils;

/**
 */
public final class TestThroughput {

  private TestThroughput() { }

  /**
   * @param args arguments
   * @throws IOException if something
   */
  public static void main(String[] args) throws IOException {
    final int t = Integer.parseInt(args[2]);
    final SimpleThreadPool stp = new SimpleThreadPool(t, "blah", true);
    for (int i = 0; i < t; ++i) {
      stp.execute(new Foo(args));
    }
    stp.terminate();
  }

  private static final class Foo implements IORunnable {
    private final String[] mArgs;
    private Foo(String[] args) {
      mArgs = args;
    }

    @Override
    public void run() throws IOException {

      final int numRecs = Integer.parseInt(mArgs[0]);
      final TempRecordWriter wrt;
      final File f = File.createTempFile("boo", "yah");
      wrt = new TempRecordWriterNio(FileUtils.createOutputStream(f, true));
      //   final Random r = new Random(75521593);
      for (int i = 0; i < numRecs; ++i) {
        final BinaryTempFileRecord rec = new BinaryTempFileRecord(true, false, false, false);
        rec.setStartPosition(i);
        rec.setReadId(i);
        rec.setReferenceId(i);
        rec.setAlignmentScore(i);
        rec.setComboScore(i);
        rec.setMatePosition(i);
        rec.setTemplateLength(i);
        rec.setNumberMismatches(i);
        rec.setSamFlags((byte) (i % 256));
        final byte[] ciggy = new byte[i % 30];
        rec.setCigarString(ciggy);
        final byte[] md5y = new byte[i % 30];
        rec.setMdString(md5y);
        wrt.writeRecord(rec);
      }
      wrt.close();
      final TempRecordReader rd;
      rd = new TempRecordReaderNio(FileUtils.createGzipInputStream(f, false), new TempRecordReader.RecordFactory(true, false, false, false));
      for (int i = 0; i < numRecs; ++i) {
        final BinaryTempFileRecord rec = rd.readRecord();
        assert !rec.isSentinelRecord();
      }
      rd.close();
      if (!f.delete()) {
        throw new IOException("Could not delete file: " + f.getPath());
      }
    }
  }
}
