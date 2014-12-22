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
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 * Test class
 */
public class FileCompressedOutputStreamTest extends TestCase {

  public void testSomeMethod() throws IOException {
    final File dir = FileUtils.createTempDir("cstreamtest", "test");
    try {
      PortableRandom pr = new PortableRandom(42);
      byte[] b = new byte[1000000];
      for (int i = 0; i < b.length; i++) {
        b[i] = (byte) pr.nextInt(CompressedMemorySequencesReader.MAX_QUAL_VALUE);
      }
      CompressedByteArray bwba = new CompressedByteArray(b.length, CompressedMemorySequencesReader.MAX_QUAL_VALUE, false);
      bwba.set(0, b, b.length);
      File f = new File(dir, "cstreamtest");
      try (DataOutputStream dos = new DataOutputStream(new FileOutputStream(f))) {
        bwba.dumpCompressedValues(dos, b.length);
      }
      File fo = new File(dir, "costreamtest");
      try (FileCompressedOutputStream out = new FileCompressedOutputStream(fo, CompressedMemorySequencesReader.MAX_QUAL_VALUE)) {
        for (int i = 0; i < b.length; ) {
          int amount = pr.nextInt(b.length - i + 1);
          out.write(b, i, amount);
          i += amount;
          if (pr.nextInt(4) < 1) {
            out.flush();
          }
        }
      }
      byte[] exp = FileBitwiseOutputStreamTest.fileToByteArray(f);
      byte[] res = FileBitwiseOutputStreamTest.fileToByteArray(fo);
      assertEquals(exp.length, res.length);
      for (int i = 0; i < exp.length; i++) {
        assertEquals("i: " + i + " exp: " + exp[i] + " res: " + res[i], exp[i], res[i]);
      }
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
    //assertTrue(Arrays.equals(exp, res));
  }

  public void testCanRead() throws IOException {
    final File dir = FileUtils.createTempDir("bwstreamtest", "test");
    try {
      PortableRandom pr = new PortableRandom(42);
      byte[] b = new byte[1000000];
      for (int i = 0; i < b.length; i++) {
        b[i] = (byte) pr.nextInt(CompressedMemorySequencesReader.MAX_QUAL_VALUE);
      }
      final File outF = new File(dir, "out");
      final FileCompressedOutputStream fos = new FileCompressedOutputStream(outF, CompressedMemorySequencesReader.MAX_QUAL_VALUE);
      try {
        fos.write(b);
      } finally {
        fos.close();
      }
      byte[] res = new byte[b.length];
      try (FileCompressedInputStream fis = new FileCompressedInputStream(outF, CompressedMemorySequencesReader.MAX_QUAL_VALUE, fos.values(), false)) {
        int pos = 0;
        int len;
        while ((len = fis.read(res, pos, res.length - pos)) > 0) {
          pos += len;
        }
      }
      assertEquals(b.length, res.length);
      for (int i = 0; i < b.length; i++) {
        assertEquals("i: " + i + " exp: " + b[i] + " res: " + res[i], b[i], res[i]);
      }
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

}
