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
package com.rtg.position.output;

import java.io.IOException;

import com.rtg.launcher.BuildParams;
import com.rtg.launcher.ISequenceParams;
import com.rtg.launcher.MockReaderParams;
import com.rtg.launcher.MockSequenceParams;
import com.rtg.launcher.ReaderParams;
import com.rtg.mode.SequenceMode;
import com.rtg.util.diagnostic.Diagnostic;



/**
 *
 */
public class GapBucketsInfoTest extends GapBucketsTest {

  public void test2BucketEx() throws IOException {
    Diagnostic.setLogStream();
    final int[] seqLengths = {0, 1, 7, 8, 9, 15, 16, 17};
    final GapBuckets<Object> bu = makeBuckets(8, 0, seqLengths, 1, SequenceMode.BIDIRECTIONAL);

    //testing isempty
    assertTrue(bu.isEmpty());
    bu.set(0, new Object());
    assertFalse(bu.isEmpty());
    bu.set(0, null);
    assertTrue(bu.isEmpty());
    bu.set(3, new Object());
    assertFalse(bu.isEmpty());

    //detailed check of bucketEx
    assertEquals(7, bu.bucketEx(3, 0, 0));
    assertEquals(7, bu.bucketEx(3, 0, 0));
    assertEquals(6, bu.bucketEx(3, 0, 1));
    assertEquals(6, bu.bucketEx(3, 0, 2));
    assertEquals(6, bu.bucketEx(3, 0, 8));
    assertEquals(5, bu.bucketEx(3, 0, 9));
    assertEquals(7, bu.bucketEx(3, 7, 0));
    assertEquals(8, bu.bucketEx(3, 8, 0));
    assertEquals(8, bu.bucketEx(3, 15, 0));
    assertEquals(9, bu.bucketEx(3, 16, 0));
    assertEquals(9, bu.bucketEx(3, 23, 0));
    assertEquals(10, bu.bucketEx(3, 24, 0));
    assertEquals(6, bu.bucketEx(3, -1, 0));

  }


  public void test3() throws IOException {
    Diagnostic.setLogStream();
    final int[] seqLengths = {0, 1, 7, 8, 9, 15, 16, 17};
    final GapBuckets<Object> bu = makeBuckets(8, 0, seqLengths, 1, SequenceMode.BIDIRECTIONAL);
    bu.globalIntegrity();

    //detailed check of bucket
    int j = 0;
    for (int i24 = 0; i24 < 4 * 24; i24 += 24) {
      for (int i8 = 0, k = 0; i8 < 24; i8 += 8, ++k) {
        for (int i = 0; i < 8; ++i) {
          assertEquals(k + 7, bu.bucket(3, j, 0));
          assertEquals(9 - k, bu.bucket(3, 0, j + 1));
          ++j;
        }
      }
    }
  }

  /**
   * Non-zero start.
   */
//  public void testNonZeroStart() {
//    Diagnostic.setLogStream();
//    final int[] seqLengths = {0, 1, 7, 8, 9, 15, 16, 17};
//    try {
//      makeBuckets(8, 1, seqLengths, 1, SequenceMode.UNIDIRECTIONAL);
//      fail();
//    } catch (final RuntimeException e) {
//      assertEquals("Unable to deal with non-zero start position.", e.getMessage());
//    }
//  }

  public void testManyBuckets() throws IOException {
    Diagnostic.setLogStream();
    final int[] seqLengths = {(Integer.MAX_VALUE >> 1) + 2, Integer.MAX_VALUE >> 1};
    makeBuckets(1, 0, seqLengths, 1, SequenceMode.UNIDIRECTIONAL);
  }

  public void testTooManySequences() throws IOException {
    Diagnostic.setLogStream();
    final SequenceMode mode = SequenceMode.UNIDIRECTIONAL;
    final ReaderParams srp = new MockReaderParams(0, Integer.MAX_VALUE + 1L, mode.codeType());
    final ISequenceParams subjectParams = new MockSequenceParams(srp, mode, 0, 0);
    final BuildParams params = BuildParams.builder().windowSize(1).stepSize(1).sequences(subjectParams).create();
    try {
      new GapBucketsInfo(params, 0, 1, 1000);
      fail();
    } catch (final RuntimeException e) {
      assertEquals("too many sequences.2147483648", e.getMessage());
    }
  }

}
