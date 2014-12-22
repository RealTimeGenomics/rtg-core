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

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

import com.rtg.util.PortableRandom;
import com.rtg.util.gzip.WorkingGzipInputStream;

import junit.framework.TestCase;

/**
 */
public class AdjustableGZIPOutputStreamTest extends TestCase {

  private class DummyOutputStream extends AdjustableGZIPOutputStream {

    public DummyOutputStream(OutputStream out) throws IOException {
      super(out);
    }

    public void checkBufferSize(int bufferSize) {
      assertEquals(bufferSize, buf.length);
    }
  }
  public void testBufSize() throws IOException {
    final ByteArrayOutputStream bos = new ByteArrayOutputStream();
    final DummyOutputStream gzipStream = new DummyOutputStream(bos);
    gzipStream.checkBufferSize(65536);
  }


  public void testCompression() throws IOException {
    // generate the uncompressed file/contents.
    final String contents = makeJunk();

    // now check that levels 1..3 of compression give decreasing sizes.
    long size = contents.length();
    long size2 = 0;
    for (int level = 1; level <= 3; level++) {
      final long compressedSize = compress(contents, level);
      //System.out.println("level " + level + " compresses to " + compressedSize + " bytes");
      assertTrue(compressedSize < size);
      size = compressedSize;
      if (level == 2) {
        size2 = compressedSize;
      }
    }

    // now check the default constructor does level 2.
    final ByteArrayOutputStream bos2 = new ByteArrayOutputStream();
    final OutputStream out = new AdjustableGZIPOutputStream(bos2);
    out.write(contents.getBytes());
    out.flush();
    out.close();
    final long size3 = bos2.size();
    assertEquals(size2, size3);
  }

  /**
   * Compress <code>contents</code> with the given level of compression.
   * @param level compression level (1..9)
   * @return the size of the compressed contents in bytes.
   * @throws IOException
   */
  public long compress(String contents, int level) throws IOException {
    final ByteArrayOutputStream bos = new ByteArrayOutputStream();
    final OutputStream out = new AdjustableGZIPOutputStream(bos, 1024, level);
    out.write(contents.getBytes());
    out.flush();
    out.close();
    final byte[] bytes = bos.toByteArray();
    final long size = bos.size();

    // now check that decompressing gives the original contents
    final InputStream in = new WorkingGzipInputStream(new ByteArrayInputStream(bytes));
    for (int i = 0; i < contents.length(); i++) {
      assertEquals(contents.charAt(i), in.read());
    }
    assertEquals(-1, in.read());
    return size;
  }

  private String makeJunk() {
    final String[] words = {"the", "theme", "of", "this", "is", "to", "have",
        "a", "medium", "level", "of", "randomness!"};
    final StringBuilder builder = new StringBuilder();
    final PortableRandom ran = new PortableRandom(42);
    for (int i = 0; i < 1000; i++) {
      final String word = words[ran.nextInt(words.length)];
      builder.append(word);
    }
    return builder.toString();
  }
}
