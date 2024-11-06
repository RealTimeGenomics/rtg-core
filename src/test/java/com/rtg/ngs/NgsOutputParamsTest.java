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
package com.rtg.ngs;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.OutputStream;

import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 */
public class NgsOutputParamsTest extends TestCase {

  @Override
  protected void setUp() throws Exception {
    super.setUp();
    Diagnostic.setLogStream();
  }

  NgsOutputParams getParams(final int topN, final int errorLimit, final boolean progress, final String logFile, final boolean exclude, final boolean useIds, final boolean tabular, final boolean sorted, final int numFiles) {
    final NgsFilterParams filterParams = NgsFilterParams.builder().outputFilter(OutputFilter.NONE).topN(topN).exclude(exclude).useids(useIds).errorLimit(errorLimit).create();
    return NgsOutputParams.builder().progress(progress).outputDir(new File(logFile)).filterParams(filterParams).tabular(tabular).sorted(sorted).numberOfOutputFiles(numFiles).create();
  }

  public void testEquals() {
    final NgsOutputParams a1 = getParams(10, 5, false, "log", false, false, true, true, 1);
    final NgsOutputParams a2 = getParams(10, 5, false, "log", false, false, true, true, 1);
    final NgsOutputParams c = getParams(10, 5, false, "foo", false, false, true, true, 1);
    final NgsOutputParams d = getParams(10, 3, false, "log", false, false, true, true, 1);
    final NgsOutputParams e = getParams(10, 5, true, "log", false, false, true, true, 1);
    final NgsOutputParams f = getParams(11, 5, true, "log", false, false, true, true, 1);
    final NgsOutputParams g = getParams(10, 5, true, "log", true, false, true, true, 1);
    final NgsOutputParams h = getParams(10, 5, true, "log", true, true, true, true, 1);
    assertTrue(h.tabular());
    TestUtils.equalsHashTest(new NgsOutputParams[][] {{a1, a2}, {c}, {d}, {e}, {f}, {g}, {h}});
  }

  public void testStreams() throws IOException {
    final File dir = File.createTempFile("NgsOutputParams", "streams");
    assertTrue(dir.delete());
    final NgsOutputParams params = getParams(10, 5, true, dir.getPath(), false, false, true, true, 1);
    try {
      final OutputStream repeat = params.repeatsStream();
      checkStream(repeat, new File(dir, "ambiguous"));
    } finally {
      assertTrue(!dir.exists() || FileHelper.deleteAll(dir));
    }
    try {
      final OutputStream unmapped = params.unmappedStream();
      checkStream(unmapped, new File(dir, "unmapped"));
    } finally {
      assertTrue(!dir.exists() || FileHelper.deleteAll(dir));
    }
  }

  private void checkStream(final OutputStream stream, final File file) throws IOException {
    final byte[] testArr = {0, 1, 2, 3, 4, 5};
    try {
      stream.write(testArr);
    } finally {
      stream.close();
    }
    try (FileInputStream fis = new FileInputStream(file)) {
      for (final byte a : testArr) {
        assertEquals((int) a, fis.read());
      }
      assertEquals(-1, fis.read());
    }
  }

  public void test0() {
    final NgsOutputParams sp = getParams(10, 5, true, "Foo", false, false, false, true, 1);
    sp.integrity();
    assertEquals(5, sp.errorLimit());
    assertEquals(true, sp.progress());
    assertEquals(new File(new File("Foo"), "zorb"), sp.file("zorb"));
    assertEquals(OutputFilter.NONE, sp.outFilter());
    assertEquals(10, sp.topN());
    assertEquals(false, sp.exclude());
    assertEquals(""
        + "NgsOutputParams progress=" + Boolean.TRUE.toString()
        + " output=Foo filterParams={NgsFilterParams filter=NONE topN=10 maxTopResults=5 error limit=5 mated max mismatches=10% unmated max mismatches=10% exclude=" + Boolean.FALSE.toString()
        + " use-ids=" + Boolean.FALSE.toString()
        + " zip=" + Boolean.FALSE.toString()
        + "} tabular=" + Boolean.FALSE.toString()
        + " sorted=" + Boolean.TRUE.toString() + " tempDir=" + null + " numoutfiles=1 bam=" + Boolean.FALSE.toString() + " sam=" + Boolean.TRUE.toString() + " sdf=" + Boolean.FALSE.toString() + " keepIntermediate=" + Boolean.FALSE.toString() + " mergeMatchResults=" + Boolean.TRUE.toString() + " mergeAlignmentResults=" + Boolean.TRUE.toString()
        , sp.toString()
    );

    assertFalse(sp.isCompressOutput());
    assertEquals(1, sp.numberOfOutputFiles());

    final NgsOutputParams sp2 = getParams(10, 5, true, "Foo", false, false, false, true, 123);
    assertEquals(123, sp2.numberOfOutputFiles());
  }

  public void testMisc() {
    final NgsOutputParamsBuilder builder = new NgsOutputParamsBuilder();
    assertFalse(builder.mBam);
    assertFalse(builder.mProgress);
    assertFalse(builder.mTabular);
    assertFalse(builder.mSorted);
    assertEquals(1, builder.mNumOutputFiles);
    assertFalse(builder.mKeepIntermediate);
    assertEquals(builder, builder.tempFilesDir(null));
    assertNull(builder.mTempFilesDir);
    assertEquals(builder, builder.bam(true));
    assertEquals(builder, builder.keepIntermediateFiles(true));
  }

}

