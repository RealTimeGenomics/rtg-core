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

import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;

import com.rtg.util.PortableRandom;
import com.rtg.util.bytecompression.CompressedByteArray;

import junit.framework.TestCase;

/**
 * Test class
 */
public class FileCompressedInputStreamTest extends TestCase {

  public void testSomeMethod() throws IOException {
    PortableRandom pr = new PortableRandom(42);
    byte[] b = new byte[1000000];
    for (int i = 0; i < b.length; i++) {
      b[i] = (byte) pr.nextInt(5);
    }
    CompressedByteArray bwba = new CompressedByteArray(b.length, 64, false);
    bwba.set(0, b, b.length);
    File f = File.createTempFile("bwstreamtest", "test");
    try {
      try (DataOutputStream dos = new DataOutputStream(new FileOutputStream(f))) {
        bwba.dumpCompressedValues(dos, b.length);
      }
      try (FileCompressedInputStream fbs = new FileCompressedInputStream(f, 64, b.length, true)) {
        for (int i = 0; i < 1000; i++) {
          int pos = pr.nextInt(b.length);
          fbs.seek(pos);
          assertEquals(bwba.get(pos), (byte) fbs.read());
        }
        fbs.seek(0);
        for (int i = 0; i < b.length; i++) {
          assertEquals(bwba.get(i), (byte) fbs.read());
        }
        assertEquals(-1, fbs.read());
        fbs.seek(0);
        byte[] buf = new byte[1000];
        int pos1 = 0;
        int len;
        while ((len = fbs.read(buf, 0, buf.length)) > 0) {
          for (int i = 0; i < len; i++) {
            assertEquals(bwba.get(pos1 + i), buf[i]);
          }
          pos1 += len;
        }
        assertEquals(b.length, pos1);
      }
    } finally {
      assertTrue(!f.exists() || f.delete());
    }
  }

}
