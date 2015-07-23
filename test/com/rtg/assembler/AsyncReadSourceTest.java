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

package com.rtg.assembler;

import static com.rtg.util.StringUtils.LS;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import com.rtg.mode.DnaUtils;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.util.SimpleThreadPool;
import com.rtg.util.diagnostic.Diagnostic;

import junit.framework.TestCase;

/**
 */
public class AsyncReadSourceTest extends TestCase {

  @Override
  public void setUp() {
    Diagnostic.setLogStream();
  }

  static String fragementToString(List<byte[]> fragments) {
    final StringBuilder sb = new StringBuilder();
    sb.append("[");
    String join = "";
    for (byte[] b : fragments) {
      sb.append(join);
      sb.append(DnaUtils.bytesToSequenceIncCG(b));
      join = ", ";
    }
    sb.append("]");
    return sb.toString();
  }

  static void listEquals(List<byte[]> a, List<byte[]> b) {
    assertEquals("expected <" + fragementToString(a) + "> but was <" + fragementToString(b) + ">", a.size(), b.size());
    for (int i = 0; i < a.size(); i++) {
      assertTrue("expected <" + fragementToString(a) + "> but was <" + fragementToString(b) + ">", Arrays.equals(a.get(i), b.get(i)));
    }
  }

  public void testReadPairSource() throws IOException {
    final SequencesReader left = ReaderTestUtils.getReaderDnaMemory(ReadPairSourceTest.LEFT_SEQUENCE);
    final SequencesReader middle = ReaderTestUtils.getReaderDnaMemory(">a" + LS + "AAAA" + LS + ">b" + LS + "ATAT" + LS);
    final SequencesReader right = ReaderTestUtils.getReaderDnaMemory(ReadPairSourceTest.RIGHT_SEQUENCE);
    final ReadPairSource source = new ReadPairSource(left, middle, right);
    assertEquals(2, source.numberReads());
    AsyncReadSource async = new AsyncReadSource(source, "testAsyncReadSource");
    SimpleThreadPool pool = new SimpleThreadPool(2, "foo", true);
    pool.execute(async);
    List<byte[]> fragments = async.nextFragments();
    assertNotNull(fragments);
    listEquals(Arrays.asList(DnaUtils.encodeString("ACGT"), DnaUtils.encodeString("AAAA"), DnaUtils.encodeString("CCCCAA")), fragments);

    fragments = async.nextFragments();
    assertNotNull(fragments);
    listEquals(Arrays.asList(DnaUtils.encodeString("GGGG"), DnaUtils.encodeString("ATAT"), DnaUtils.encodeString("ACAAAC")), fragments);

    assertNull(async.nextFragments());
    assertNull(async.nextFragments());

    source.reset();
    pool.terminate();
    pool = new SimpleThreadPool(2, "foo", true);
    async = new AsyncReadSource(source, "testAsyncReadSource");
    pool.execute(async);
    assertNotNull(async.nextFragments());
    fragments = async.nextFragments();
    assertNotNull(fragments);
    listEquals(Arrays.asList(DnaUtils.encodeString("GGGG"), DnaUtils.encodeString("ATAT"), DnaUtils.encodeString("ACAAAC")), fragments);

    source.close();
    pool.terminate();
  }
}
