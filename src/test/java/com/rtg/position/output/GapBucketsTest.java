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
package com.rtg.position.output;

import static com.rtg.util.StringUtils.LS;

import java.io.IOException;

import com.rtg.launcher.BuildParams;
import com.rtg.launcher.ISequenceParams;
import com.rtg.launcher.MockReaderParams;
import com.rtg.launcher.MockSequenceParams;
import com.rtg.launcher.ReaderParams;
import com.rtg.mode.SequenceMode;
import com.rtg.reader.MockArraySequencesReader;
import com.rtg.reader.SequencesReader;
import com.rtg.util.diagnostic.Diagnostic;

import junit.framework.TestCase;

/**
 */
public class GapBucketsTest extends TestCase {

  static BuildParams makeParams(final int stepSize, final int[] seqLengths, final int start, final SequenceMode mode) throws IOException {
    final SequencesReader reader = new MockArraySequencesReader(mode.type(), seqLengths);
    final ReaderParams srp = new MockReaderParams(reader);
    final ISequenceParams subjectParams = new MockSequenceParams(srp, mode, start, reader.numberSequences());
    return BuildParams.builder().windowSize(stepSize).stepSize(stepSize).sequences(subjectParams).create();
  }

  protected GapBuckets<Object> makeBuckets(final int stepSize, final int start, final int[] seqLengths, final int scale, final SequenceMode mode) throws IOException {
    final BuildParams buildParams = GapBucketsTest.makeParams(stepSize, seqLengths, start, mode);
    final GapBucketsInfo bucketInfo = new GapBucketsInfo(buildParams, scale, 16, 1000);
    return new GapBuckets<>(bucketInfo);
  }

  public void testEmpty() throws IOException {
    Diagnostic.setLogStream();
    final GapBuckets<Object> bu = makeBuckets(8, 0, new int[] {}, 1, SequenceMode.UNIDIRECTIONAL);
    bu.globalIntegrity();
    assertEquals(0, bu.numberSeqIds());
    assertEquals("GapBuckets pointers[0]" + LS + "buckets [1]" + LS, bu.toString());
    assertTrue(bu.isEmpty());
  }

  /**
   * Hammers the bucket call.
   * Makes sure a wide range of values all map to different buckets.
   */
  private void checkBuckets(final int step, final GapBuckets<Object> bu, final int[] seqLengths, final int nFrames, final int gap, final int inc) {
    for (int a = 0; a < step; ++a) {
      for (int s = 0; s < bu.numberSeqIds(); ++s) {
        for (int j = inc; j <= inc + gap; ++j) {
          long last = Long.MIN_VALUE;
          for (int i = a; i < a + seqLengths[s / nFrames]; i += step) {
            final long bucket = bu.bucket(s, i, j);
            final Object obj = new Object();
            assertEquals(null, bu.get(bucket));
            bu.set(bucket, obj);
            assertTrue(obj == bu.get(bucket));
            bu.set(bucket, null);
            assertEquals(null, bu.get(bucket));
            assertTrue("last=" + last + " bucket=" + bucket + " a=" + a + " s=" + s + " j=" + j + " i=" + i + " inc=" + inc, last != bucket);
            last = bucket;
          }
        }
      }
    }

    for (int s = 0; s < bu.numberSeqIds(); ++s) {
      assertEquals(bu.firstBucket(s), bu.bucket(s, 0, 0));
      assertEquals(bu.lastBucket(s), bu.bucket(s, 0, 1));
    }
  }

  public void test1() throws IOException {
    Diagnostic.setLogStream();
    final int[] seqLengths = {0, 1, 7, 8, 9, 15, 16, 17};
    final GapBuckets<Object> bu = makeBuckets(8, 0, seqLengths, 1, SequenceMode.UNIDIRECTIONAL);
    bu.globalIntegrity();
    assertEquals(8, bu.numberSeqIds());
    final String expected = ""
      + "GapBuckets pointers[8]" + LS
      + "[0]->0..1" + LS
      + "[1]->2..4" + LS
      + "[2]->5..7" + LS
      + "[3]->8..10" + LS
      + "[4]->11..14" + LS
      + "[5]->15..18" + LS
      + "[6]->19..22" + LS
      + "[7]->23..27" + LS
      + "buckets [32]" + LS
      ;
    assertEquals(expected, bu.toString());
    assertTrue(bu.isEmpty());
    checkBuckets(8, bu, seqLengths, 1, 16, 0);
  }

  //set scale to 2
  public void test1scale() throws IOException {
    Diagnostic.setLogStream();
    final int[] seqLengths = {0, 1, 7, 8, 9, 15, 16, 17};
    final GapBuckets<Object> bu = makeBuckets(8, 0, seqLengths, 2, SequenceMode.UNIDIRECTIONAL);
    bu.globalIntegrity();
    assertEquals(8, bu.numberSeqIds());
    final String expected = ""
      + "GapBuckets pointers[8]" + LS
      + "[0]->0..0" + LS
      + "[1]->1..2" + LS
      + "[2]->3..4" + LS
      + "[3]->5..6" + LS
      + "[4]->7..8" + LS
      + "[5]->9..10" + LS
      + "[6]->11..12" + LS
      + "[7]->13..15" + LS
      + "buckets [16]" + LS
      ;
    assertEquals(expected, bu.toString());
    assertTrue(bu.isEmpty());
    checkBuckets(8, bu, seqLengths, 1, 16, 0);
  }

  public void test2() throws IOException {
    Diagnostic.setLogStream();
    final int[] seqLengths = {0, 1, 7, 8, 9, 15, 16, 17};
    final GapBuckets<Object> bu = makeBuckets(8, 0, seqLengths, 1, SequenceMode.BIDIRECTIONAL);
    bu.globalIntegrity();
    assertEquals(16, bu.numberSeqIds());
    assertEquals(64, bu.numberBuckets());
    final String expected = ""
      + "GapBuckets pointers[16]" + LS
      + "[0]->0..1" + LS
      + "[1]->2..3" + LS
      + "[2]->4..6" + LS
      + "[3]->7..9" + LS
      + "[4]->10..12" + LS
      + "[5]->13..15" + LS
      + "[6]->16..18" + LS
      + "[7]->19..21" + LS
      + "[8]->22..25" + LS
      + "[9]->26..29" + LS
      + "[10]->30..33" + LS
      + "[11]->34..37" + LS
      + "[12]->38..41" + LS
      + "[13]->42..45" + LS
      + "[14]->46..50" + LS
      + "[15]->51..55" + LS
      + "buckets [64]" + LS
      ;
    assertEquals(expected, bu.toString());
    checkBuckets(8, bu, seqLengths, 2, 16, 0);

    //testing isempty
    assertTrue(bu.isEmpty());
    bu.set(0, new Object());
    assertFalse(bu.isEmpty());
    bu.set(0, null);
    assertTrue(bu.isEmpty());
    bu.set(3, new Object());
    assertFalse(bu.isEmpty());

    //detailed check of bucket
    assertEquals(7, bu.bucket(3, 0, 0));
    assertEquals(7, bu.bucket(3, 0, 0));
    assertEquals(9, bu.bucket(3, 0, 1));
    assertEquals(9, bu.bucket(3, 0, 8));
    assertEquals(8, bu.bucket(3, 0, 9));
    assertEquals(7, bu.bucket(3, 7, 0));
    assertEquals(8, bu.bucket(3, 8, 0));
    assertEquals(8, bu.bucket(3, 15, 0));
    assertEquals(9, bu.bucket(3, 16, 0));
    assertEquals(9, bu.bucket(3, 23, 0));
    assertEquals(7, bu.bucket(3, 24, 0));

  }

  public void test2scale() throws IOException {
    Diagnostic.setLogStream();
    final int[] seqLengths = {0, 1, 7, 8, 9, 15, 16, 17};
    final GapBuckets<Object> bu = makeBuckets(8, 0, seqLengths, 2, SequenceMode.BIDIRECTIONAL);
    bu.globalIntegrity();
    assertEquals(16, bu.numberSeqIds());
    assertEquals(32, bu.numberBuckets());
    final String expected = ""
      + "GapBuckets pointers[16]" + LS
      + "[0]->0..0" + LS
      + "[1]->1..1" + LS
      + "[2]->2..3" + LS
      + "[3]->4..5" + LS
      + "[4]->6..7" + LS
      + "[5]->8..9" + LS
      + "[6]->10..11" + LS
      + "[7]->12..13" + LS
      + "[8]->14..15" + LS
      + "[9]->16..17" + LS
      + "[10]->18..19" + LS
      + "[11]->20..21" + LS
      + "[12]->22..23" + LS
      + "[13]->24..25" + LS
      + "[14]->26..28" + LS
      + "[15]->29..31" + LS
      + "buckets [32]" + LS
      ;
    assertEquals(expected, bu.toString());
    checkBuckets(8, bu, seqLengths, 2, 16, 0);

    //testing isempty
    assertTrue(bu.isEmpty());
    bu.set(0, new Object());
    assertFalse(bu.isEmpty());
    bu.set(0, null);
    assertTrue(bu.isEmpty());
    bu.set(3, new Object());
    assertFalse(bu.isEmpty());

    //detailed check of bucket
    assertEquals(4, bu.bucket(3, 0, 0));
    assertEquals(4, bu.bucket(3, 0, 0));
    assertEquals(5, bu.bucket(3, 0, 1));
    assertEquals(5, bu.bucket(3, 0, 8));
    assertEquals(4, bu.bucket(3, 0, 9));
    assertEquals(4, bu.bucket(3, 7, 0));
    assertEquals(5, bu.bucket(3, 8, 0));
    assertEquals(5, bu.bucket(3, 15, 0));
    assertEquals(4, bu.bucket(3, 16, 0));
    assertEquals(4, bu.bucket(3, 23, 0));
    assertEquals(5, bu.bucket(3, 24, 0));
  }
}
