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
package com.rtg.index.params;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

import com.rtg.util.Constants;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 */
public class CountParamsTest extends TestCase {

  CountParams getParams(final File dir, final int topN, final int minHits, final long maxFileSize) {
    return new CountParams(dir, topN, minHits, maxFileSize, false);
  }

  CountParams getParams(final File dir, final int topN, final int minHits) {
    return new CountParams(dir, topN, minHits, false);
  }

  public void testOutputStream() throws IOException {
    Diagnostic.setLogStream();
    final File temp = FileUtils.createTempDir("temp", "countparams");
    try {
      CountParams a = new CountParams(new File(temp, "blah"), 1, 1, true);
      final File output = new File(temp, "blah" + StringUtils.FS + "blah.gz");
      OutputStream stream = a.outStream("blah");
      try {
        stream.write(256);
      } finally {
        stream.close();
      }
      final long initialSize = output.length();
      assertTrue(output.isFile());
      assertTrue(output.length() > 0);
      stream = a.outStream("blah");
      stream.close();
      assertTrue(output.isFile());
      assertTrue(output.length() < initialSize);
      a = new CountParams(output, 1, 1, false);
      try {
        a.outStream("boo");
        fail();
      } catch (final IOException e) {
        assertEquals("Unable to create directory \"" + output.getPath() + "\"", e.getMessage());
      }
    } finally {
      FileHelper.deleteAll(temp);
    }
  }

  public void testEquals() {
    final File dira = new File("a");
    final File dirb = new File("b");
    final CountParams a1 = getParams(dira, 5, 2);
    final CountParams a2 = getParams(dira, 5, 2);
    final CountParams b = getParams(dirb, 5, 2);
    final CountParams c = getParams(dira, 5, 2, Constants.MINIMUM_FILE_CHUNK_SIZE);
    final CountParams d = getParams(dirb, 3, 2);
    final CountParams e = getParams(dirb, 5, 1);
    TestUtils.equalsHashTest(new CountParams[][] {{a1, a2}, {b}, {c}, {d}, {e}});
  }

  public void test0() {
    final CountParams sp = getParams(new File("Bar"), 5, 2);
    sp.integrity();
    assertEquals("Bar", sp.directory().toString());
    assertEquals(5, sp.topN());
    assertEquals(2, sp.minHits());
    assertEquals(1000000000, sp.maxFileSize());
    assertEquals(" CountParams directory=Bar topN=5 minHits=2 max. file size=1000000000", sp.toString());
  }

  public void test() {
    final CountParams sp = getParams(new File("Bar"), 5, 2, 1001);
    sp.integrity();
    assertEquals("Bar", sp.directory().toString());
    assertEquals(1001, sp.maxFileSize());
    assertEquals(" CountParams directory=Bar topN=5 minHits=2 max. file size=1001", sp.toString());
  }

}

