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
