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

package com.rtg.reader;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

import com.rtg.util.test.RandomByteGenerator;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 * Abstract class for testing the compressed and bitwise file streams.
 */
public abstract class AbstractFileStreamRegression extends TestCase {

  private static final int RANGE = 128;
  private static final long NUM_ELEMENTS = 10L * Integer.MAX_VALUE + 9000L;

  protected abstract long calcLength(int range, long elements);
  protected abstract OutputStream createOutputStream(File file, int range) throws IOException;
  protected abstract InputStream createInputStream(File file, int range, long elements, boolean seekable) throws IOException;

  /**
   * Test the input and output streams
   * @throws IOException if an error occurs
   */
  public void testStreams() throws IOException {
    Diagnostic.setLogStream();
    doTest(RANGE, NUM_ELEMENTS);
  }

  private void doTest(int range, long elements) throws IOException {
    final File outDir = FileHelper.createTempDirectory();
    System.err.println(outDir.getPath());
    try {
      final File temp = new File(outDir, "temp.bin");
      final RandomByteGenerator value = new RandomByteGenerator(range);
      final byte[] buffer = new byte[1024];
      try (OutputStream out = createOutputStream(temp, range)) {
        for (long l = 0; l < elements; ) {
          int i = 0;
          for (; i < buffer.length && l < elements; i++, l++) {
            buffer[i] = value.nextValue();
          }
          out.write(buffer, 0, i);
        }
        out.flush();
      }
      assertTrue(temp.exists());
      assertEquals(calcLength(range, elements), temp.length());

      value.reset();
      //      final long labelChunkSize = NUM_ELEMENTS / 1000L;
//      long nextLabelOutput = labelChunkSize;
      try (InputStream in = createInputStream(temp, range, elements, false)) {
        for (long l = 0; l < elements; ) {
          final int read = in.read(buffer);
          for (int i = 0; i < read; i++, l++) {
            assertEquals(value.nextValue(), buffer[i]);
          }
//          if (l >= nextLabelOutput) {
//            System.err.println("Elements Read: " + l);
//            nextLabelOutput += labelChunkSize;
//          }
        }
      }
    } finally {
      FileHelper.deleteAll(outDir);
    }
  }
}
