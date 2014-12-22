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

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.zip.GZIPOutputStream;

import com.rtg.util.test.FileHelper;


/**
 */
public class GzipAsynchInputStreamTest extends AsynchInputStreamTest {


  @Override
  AsynchInputStream getStream(File file, String text) throws IOException {
    if (text != null) {
      FileHelper.stringToGzFile(text, file);
    }
    return new GzipAsynchInputStream(file);
  }

  public void testNullFile() throws IOException {
    try {
      new GzipAsynchInputStream((File) null);
      fail("null file exception expected");
    } catch (final IllegalArgumentException e) {
      assertEquals("File cannot be null", e.getMessage());
    }
  }

  public void testReadLarge() throws IOException {
    // after 600 experiments, these seem to be the fastest settings.
    final boolean asynch = true;
    final int bufSize = 1024 * 1024;
    final int gzipSize = 64 * 1024;
    final int readSize = 1024;
    // for (int gzipSize = 1024; gzipSize <= 64 * 1024; gzipSize *= 2) {
    //    final long t0 = System.nanoTime();
    final int kbytes = 2 * 1024;  // size of uncompressed input file in Kbytes.
    final File temp = makeLargeGzipFile(kbytes);
    //File temp = new File("/home2/marku/taskFilterSam/Test10Gb.gz");
    final File tempOut = File.createTempFile("test", "asynchOut");
    try {
      long count = 0;
      //    final long t1 = System.nanoTime();
      try (GZIPOutputStream output = new GZIPOutputStream(new FileOutputStream(tempOut))) {
        InputStream input;
        if (asynch) {
          final GzipAsynchInputStream asynchInput = new GzipAsynchInputStream(temp, bufSize - 1, gzipSize - 1);
          assertEquals(bufSize - 1, asynchInput.mQueue.maxSize());
          input = asynchInput;
        } else {
          input = FileUtils.createGzipInputStream(temp, false);
        }
        try {
          final byte[] buf = new byte[readSize];
          for (; ; ) {
            final int size = input.read(buf);
            if (size <= 0) {
              break;
            }
            //      // consume some time, to emulate processing the data
            //      byte[] dummy = {'A', 'B', 'C', 'D'};
            //      for (long i = 0; i < size; i++) {
            //        dummy[(int) i & 0x3] ^= 0x01;
            //      }
            output.write(buf, 0, size);
            count += size;
          }
        } finally {
          input.close();
        }
      }
      //    final long t2 = System.nanoTime();
      assertEquals((long) kbytes * 1024, count);
      //    System.out.println("createtime=" + (t1 - t0) / 1000000.0 + ", "
      //        + gzipSize + ", " + bufSize + ", " + readSize + ", "
      //        + (asynch ? "a" : "") + "synch, " + (t2 - t1) / 1000000.0);
      //}
    } finally {
      assertTrue(FileHelper.deleteAll(tempOut));
      assertTrue(FileHelper.deleteAll(temp));
    }
  }

  /**
   * Creates a largish gzipped file containing <code>numCopies</code>
   * of <code>line</code>.
   *
   * @param kbytes size of file
   * @return the temporary File that has been created.
   * @throws IOException
   */
  public File makeLargeGzipFile(int kbytes) throws IOException {
    final byte[] line = new byte[1024];
    final File file = File.createTempFile("test", "GzipAsynch.gz");
    // fill file with a larger amount of data.
    try (GZIPOutputStream out = new GZIPOutputStream(new FileOutputStream(file))) {
      for (int i = 0; i < kbytes; i++) {
        // generate a randomish line.
        for (int pos = 0; pos < 1024; pos++) {
          line[pos] = (byte) ('!' + (pos * (long) i) % 91);
        }
        line[1023] = (byte) '\n';
        out.write(line);
      }
    }
    return file;
  }

  @Override
  public void testEmpty() {
    try {
      super.testEmpty();
      fail();
    } catch (IOException e) {
      // We require that empty output files contain 0 bytes
      // (i.e. not even including empty gzip blocks) for correct
      // fast concatenation of multi-part output files.

      // Java IOException in this case, .NET does not.

      // We will normally avoid reading such files by checking file length first.
    }
  }

  /*
  public void testConstantly() {
    byte[] fiftyMeg = new byte[500 * 1024 * 1024];
    final File tempFile = File.createTempFile("cont", "ziptest");
    try {
      final PortableRandom rand = new PortableRandom();
      rand.nextBytes(fiftyMeg);
      final GZIPOutputStream gzipOut = new GZIPOutputStream(new FileOutputStream(tempFile));
      try {
        gzipOut.write(fiftyMeg, 0, fiftyMeg.length);
      } finally {
        gzipOut.close();
      }
      final byte[] buf = new byte[2* 1024 * 1024];
      for (int i = 0; i < 1000; i++) {
        GzipAsynchInputStream asyncIn = new GzipAsynchInputStream(tempFile);
        final long start = System.currentTimeMillis();
        try {
          int len;
          while ((len = asyncIn.read(buf, 0, buf.length)) > 0) {


          }
        } finally {
          asyncIn.close();
        }
        System.out.println("Time: " + (System.currentTimeMillis() - start) + "ms");
      }
    } finally {
      assertTrue(tempFile.delete());
    }
  }
   */
}
