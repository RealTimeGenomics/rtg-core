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
package com.rtg.launcher;

import static com.rtg.util.StringUtils.LS;
import static com.rtg.util.StringUtils.TAB;

import java.io.IOException;
import java.io.Reader;
import java.io.StringReader;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.rtg.mode.SequenceType;
import com.rtg.reader.MockSequencesReader;
import com.rtg.reader.ReaderUtils;
import com.rtg.reference.ReferenceGenome;
import com.rtg.reference.Sex;

import junit.framework.TestCase;

/**
 */
public class HashingRegionTest extends TestCase {


  public void testNone() {
    assertTrue(HashingRegion.NONE.isInRange(Long.MAX_VALUE, Long.MAX_VALUE));
    assertTrue(HashingRegion.NONE.isInRange(Long.MIN_VALUE, Long.MIN_VALUE));
    assertTrue(HashingRegion.NONE.isInRange(1, 1));
    assertTrue(HashingRegion.NONE.isInRange(0, 0));
    assertTrue(HashingRegion.NONE.isInRange(-1, -1));
  }

  public void testClipping() {
    final HashingRegion region = new HashingRegion(43, 32, 90, 44, -1, -1);

    assertEquals(43, region.getStart());
    assertEquals(32, region.getStartClipPosition());
    assertEquals(90, region.getEnd());
    assertEquals(44, region.getEndClipPosition());
    assertEquals(91, region.getExclusiveEndId());


    assertEquals("[(43:32), (90:44)]", region.toString());

    // earlier template
    assertFalse(region.isInRange(42, 30));
    assertFalse(region.isInRange(42, 36));
    assertFalse(region.isInRange(42, 49));

    // start template
    assertFalse(region.isInRange(43, 31));
    assertTrue(region.isInRange(43, 32));
    assertTrue(region.isInRange(43, 39));

    // template in the middle
    assertTrue(region.isInRange(44, 30));
    assertTrue(region.isInRange(44, 32));
    assertTrue(region.isInRange(44, 40));
    assertTrue(region.isInRange(44, 44));
    assertTrue(region.isInRange(44, 55));


    // end Template
    assertTrue(region.isInRange(90, 43));
    assertFalse(region.isInRange(90, 44));
    assertFalse(region.isInRange(90, 46));

    // later template
    assertFalse(region.isInRange(91, 30));
    assertFalse(region.isInRange(91, 36));
    assertFalse(region.isInRange(91, 49));


    assertTrue(new HashingRegion(-1, -1, -1, -1, -1, -1).isInRange(-1, -1));

    final HashingRegion templateRegion = new HashingRegion(10, 20);

    assertFalse(templateRegion.isInRange(9, 1));
    assertTrue(templateRegion.isInRange(10, 1));
    assertTrue(templateRegion.isInRange(19, 1));
    assertFalse(templateRegion.isInRange(20, 1));
    assertFalse(templateRegion.isInRange(21, 1));
    assertEquals(20, templateRegion.getExclusiveEndId());


    // Test NONE
    assertTrue(HashingRegion.NONE.isInRange(-1, -1));
    assertTrue(HashingRegion.NONE.isInRange(10));
    assertEquals("[All inclusive]", HashingRegion.NONE.toString());



  }

  public void testTemplateClipping() {
    HashingRegion region = new HashingRegion(3, 5);

    assertFalse(region.isInRange(2));
    assertTrue(region.isInRange(3));
    assertTrue(region.isInRange(4));
    assertFalse(region.isInRange(5));

    region = new HashingRegion(3, 3, 5, 5, -1, -1);

    assertFalse(region.isInRange(2));
    assertTrue(region.isInRange(3));
    assertTrue(region.isInRange(4));
    assertTrue(region.isInRange(5));
    assertFalse(region.isInRange(6));
  }
  HashingRegion[] splitDefault(int[] sequenceLengths, int start, int numberChunks, int padding) {
    return HashingRegion.splitWorkload(sequenceLengths, start, numberChunks, HashingRegion.DEFAULT_MIN_CHUNK_SIZE, padding);
  }

  public void testSplitWorkload() {
    final HashingRegion[] ranges = splitDefault(new int[]{1000000, 1000000, 1000000, 1000000}, 0, 4, 10);

    assertEquals(4, ranges.length);
    assertEquals(0, ranges[0].getStart());
    assertEquals(0, ranges[0].getStartClipPosition());
    assertEquals(0, ranges[0].getEnd());
    assertEquals(1000000, ranges[0].getEndClipPosition());

    assertEquals(0, ranges[1].getStart());
    assertEquals(1000000, ranges[1].getStartClipPosition());
    assertEquals(1, ranges[1].getEnd());
    assertEquals(1000000, ranges[1].getEndClipPosition());

    assertEquals(1, ranges[2].getStart());
    assertEquals(1000000, ranges[2].getStartClipPosition());
    assertEquals(2, ranges[2].getEnd());
    assertEquals(1000000, ranges[2].getEndClipPosition());

    assertEquals(2, ranges[3].getStart());
    assertEquals(1000000, ranges[3].getStartClipPosition());
    assertEquals(3, ranges[3].getEnd());
    assertEquals(1000000, ranges[3].getEndClipPosition());
  }

  public void testSplitWorkload2() {
    final HashingRegion[] ranges = splitDefault(new int[]{500000, 500000, 2000000, 1000000}, 0, 4, 10);

    assertEquals(4, ranges.length);
    assertEquals(0, ranges[0].getStart());
    assertEquals(0, ranges[0].getStartClipPosition());
    assertEquals(1, ranges[0].getEnd());
    assertEquals(500000, ranges[0].getEndClipPosition());

    assertEquals(1, ranges[1].getStart());
    assertEquals(500000, ranges[1].getStartClipPosition());
    assertEquals(2, ranges[1].getEnd());
    assertEquals(1000000, ranges[1].getEndClipPosition());

    assertEquals(2, ranges[2].getStart());
    assertEquals(1000000, ranges[2].getStartClipPosition());
    assertEquals(2, ranges[2].getEnd());
    assertEquals(2000000, ranges[2].getEndClipPosition());

    assertEquals(2, ranges[3].getStart());
    assertEquals(2000000, ranges[3].getStartClipPosition());
    assertEquals(3, ranges[3].getEnd());
    assertEquals(1000000, ranges[3].getEndClipPosition());
    checkRangeArray(ranges);
  }

  public void testSplitWorkload3() {
    final HashingRegion[] ranges = splitDefault(new int[]{5689, 5329401, 2023432, 1323232, 34342, 2344343}, 0, 9, 10);

    assertEquals(9, ranges.length);
    assertEquals(0, ranges[0].getStart());
    assertEquals(0, ranges[0].getStartClipPosition());

    assertEquals(5, ranges[8].getEnd());
    assertEquals(2344343, ranges[8].getEndClipPosition());
    checkRangeArray(ranges);
  }

  public void testSplitWorkload4() {
    final HashingRegion[] ranges = HashingRegion.splitWorkload(new int[]{50, 100, 200, 400, 650, 150}, 0, 10, 10, 10);

    assertEquals(10, ranges.length);
    checkRange(ranges[0], 0, 0, 0, 2, 5, 15);
    checkRange(ranges[1], 2, 5, 0, 2, 160, 170);
    checkRange(ranges[2], 2, 160, 150, 3, 115, 125);
    checkRange(ranges[3], 3, 115, 105, 3, 270, 280);
    checkRange(ranges[4], 3, 270, 260, 4, 25, 35);
    checkRange(ranges[5], 4, 25, 15, 4, 180, 190);
    checkRange(ranges[6], 4, 180, 170, 4, 335, 345);
    checkRange(ranges[7], 4, 335, 325, 4, 490, 500);
    checkRange(ranges[8], 4, 490, 480, 4, 645, 650);
    checkRange(ranges[9], 4, 645, 635, 5, 150, 150);

    assertEquals(0, ranges[3].getReferenceStart(10, 10));
    assertEquals(0, ranges[1].getReferenceStart(2, 10));
    assertEquals(170, ranges[1].getReferenceEnd(2, 10, 200));
    assertEquals(105, ranges[3].getReferenceStart(3, 10));
    assertEquals(280, ranges[3].getReferenceEnd(3, 10, 400));

    checkRangeArray(ranges);
  }

  public void testSplitWorkLoadSpecial() {
    final HashingRegion[] ranges = HashingRegion.splitWorkload(new int[]{1, 3}, 0, 5, 10, 1);
    assertEquals(4, ranges.length);
    checkRange(ranges[0], 0, 0, 0, 0, 1, 1);
    checkRange(ranges[1], 1, 0, 0, 1, 1, 2);
    checkRange(ranges[2], 1, 1, 0, 1, 2, 3);
    checkRange(ranges[3], 1, 2, 1, 1, 3, 3);
  }

  public void testSplitWorkLoad5() {
    final HashingRegion[] ranges = HashingRegion.splitWorkload(new int[]{10, 12}, 0, 2, 10, 10);
    assertEquals(2, ranges.length);
    checkRange(ranges[0], 0, 0, 0, 1, 1, 11);
    checkRange(ranges[1], 1, 1, 0, 1, 12, 12);
  }

  public void testSplitWorkLoadMinChunkSize() {
    //Testing min chunk size
    final HashingRegion[] ranges = HashingRegion.splitWorkload(new int[]{15}, 0, 3, 6, 1);
    assertEquals(2, ranges.length);
    checkRange(ranges[0], 0, 0, 0, 0, 7, 8);
    checkRange(ranges[1], 0, 7, 6, 0, 15, 15);
  }

  private static void checkRange(HashingRegion r, long startId, long startPos, long startPad, long endId, long endPos, long endPad) {
    assertEquals(startId, r.getStart());
    assertEquals(startPos, r.getStartClipPosition());
    assertEquals(startPad, r.getStartPaddedPosition());
    assertEquals(endId, r.getEnd());
    assertEquals(endPos, r.getEndClipPosition());
    assertEquals(endPad, r.getEndPaddedPosition());
  }
//
//  public void testSplitWorkloadLongHandling() {
//    final long length = ((long) Integer.MAX_VALUE) + 4;
//    final SpanningRegion[] ranges = SpanningRegion.splitWorkload(new long[] {length, length}, 0, 2, 5);
//
//    assertEquals(2, ranges.length);
//
//    assertEquals(0, ranges[0].getStart());
//    assertEquals(0, ranges[0].getStartClipPosition());
//    assertEquals(0, ranges[0].getEnd());
//    assertEquals(length, ranges[0].getEndClipPosition());
//
//    assertEquals(0, ranges[1].getStart());
//    assertEquals(length, ranges[1].getStartClipPosition());
//    assertEquals(1, ranges[1].getEnd());
//    assertEquals(length, ranges[1].getEndClipPosition());
//    checkRangeArray(ranges);
//  }

  public void testSplitWorkloadSmall() {
    final HashingRegion[] ranges = splitDefault(new int[]{2, 2}, 0, 8, 5);

    assertEquals(4, ranges.length);

    assertEquals(0, ranges[0].getStart());
    assertEquals(0, ranges[0].getStartClipPosition());
    assertEquals(0, ranges[0].getEnd());
    assertEquals(1, ranges[0].getEndClipPosition());

    assertEquals(0, ranges[1].getStart());
    assertEquals(1, ranges[1].getStartClipPosition());
    assertEquals(0, ranges[1].getEnd());
    assertEquals(2, ranges[1].getEndClipPosition());

    assertEquals(1, ranges[3].getStart());
    assertEquals(1, ranges[3].getStartClipPosition());
    assertEquals(1, ranges[3].getEnd());
    assertEquals(2, ranges[3].getEndClipPosition());
  }

  private void checkRangeArray(HashingRegion[] regions) {
    for (int i = 1; i < regions.length; ++i) {
      if (regions[i] != null) {
        assertTrue("Clip Region end id should be the start id of the next region", regions[i].getStart() == regions[i - 1].getEnd());
        assertTrue("Clip Region end position should be the start position of the next region", regions[i].getStartClipPosition() == regions[i - 1].getEndClipPosition());
        assertTrue(regions[i].getEnd() >= regions[i - 1].getEnd());
        if (regions[i].getEnd() == regions[i - 1].getEnd()) {
          assertTrue("ClipRegions shouldn't have 0 length", regions[i].getEndClipPosition() >= regions[i - 1].getEndClipPosition());
        }
      }
    }
  }

  public void testEquals() {
    final HashingRegion r = new HashingRegion(1, 2, 3, 4, -1, -1);
    assertTrue(r.equals(r));
    assertTrue(r.equals(new HashingRegion(1, 2, 3, 4, -1, -1)));

    assertFalse(r.equals(new Object()));
    assertFalse(r.equals(null));

    assertFalse(r.equals(new HashingRegion(2, 2, 3, 4, -1, -1)));
    assertFalse(r.equals(new HashingRegion(1, 3, 3, 4, -1, -1)));
    assertFalse(r.equals(new HashingRegion(1, 2, 4, 4, -1, -1)));
    assertFalse(r.equals(new HashingRegion(1, 2, 3, 5, -1, -1)));

    assertFalse(new HashingRegion(-1, -1, -1, -1, -1, -1).equals(HashingRegion.NONE));

    assertTrue(new HashingRegion(5, -1, 6, -1, -1, -1).equals(new HashingRegion(5, 6)));
  }

  public void testHashCode() {
    assertEquals(36352083, new HashingRegion(0, 0, 0, 0, -1, -1).hashCode());
    final HashingRegion r = new HashingRegion(1L << 32 + 1, 1L << 32 + 1, 1L << 32 + 1, 1L << 32 + 1, -1, -1);
    assertEquals(36769923, r.hashCode());
  }


  public void testCompare() {
    final HashingRegion a = new HashingRegion(1, 2, 3, 4, -1, -1);
    final HashingRegion b = new HashingRegion(1, 2, 3, 5, -1, -1);
    final HashingRegion c = new HashingRegion(1, 2, 4, 4, -1, -1);
    final HashingRegion d = new HashingRegion(1, 3, 3, 4, -1, -1);
    final HashingRegion e = new HashingRegion(2, 2, 3, 4, -1, -1);

    assertEquals(0, a.compareTo(a));
    assertTrue(a.compareTo(null) > 0);

    assertTrue(a.compareTo(b) < 0);
    assertTrue(a.compareTo(c) < 0);
    assertTrue(a.compareTo(d) < 0);
    assertTrue(a.compareTo(e) < 0);

    assertTrue(b.compareTo(a) > 0);
    assertTrue(c.compareTo(a) > 0);
    assertTrue(d.compareTo(a) > 0);
    assertTrue(e.compareTo(a) > 0);
  }

  public void testLongestSubsequence() throws IOException {
    final FakeReader reader = new FakeReader(new int[] {100, 200, 300});
    HashingRegion r = new HashingRegion(1, 1);
    assertEquals(0, r.longestSubSequence(reader));
    r = new HashingRegion(0, 1);
    assertEquals(100, r.longestSubSequence(reader));
    r = new HashingRegion(0, 3);
    assertEquals(300, r.longestSubSequence(reader));
    r = new HashingRegion(0, 5, 0, 95, 0, 100);
    assertEquals(100, r.longestSubSequence(reader));
    r = new HashingRegion(0, 15, 1, -1, 5, 100);
    assertEquals(95, r.longestSubSequence(reader));
    r = new HashingRegion(0, 15, 0, 100, 5, 100);
    assertEquals(95, r.longestSubSequence(reader));
    r = new HashingRegion(0, 30, 0, 70, 20, 80);
    assertEquals(60, r.longestSubSequence(reader));
    r = new HashingRegion(0, 0, 2, 250, 0, 260);
    assertEquals(260, r.longestSubSequence(reader));
    r = new HashingRegion(0, 0, 2, 180, 0, 190);
    assertEquals(200, r.longestSubSequence(reader));
    r = new HashingRegion(0, 20, 1, 40, 10, 50);
    assertEquals(90, r.longestSubSequence(reader));
    assertEquals(300, HashingRegion.NONE.longestSubSequence(reader));
  }

  private static final class FakeReader extends MockSequencesReader {
    private final int[] mLengths;
    FakeReader(int[] lengths) {
      super(SequenceType.DNA);
      mLengths = lengths;
    }

    @Override
    public int[] sequenceLengths(long start, long end) {
      final int[] ret = new int[(int) (end - start)];
      System.arraycopy(mLengths, (int) start, ret, 0, ret.length);
      return ret;
    }

    @Override
    public long maxLength() {
      long max = 0;
      for (int mLength : mLengths) {
        if (max < mLength) {
          max = mLength;
        }
      }
      return max;
    }
  }

  public void testNewSplitWorkload() throws IOException {
    final MockSequencesReader r = new MockSequencesReader(SequenceType.DNA, 3);
    r.setLengths(new int[] {10000000, 20000000, 30000000});
    final String expected = ""
            + "[[(0:0), (0:3000000)]"
            + ", [(0:3000000), (0:6000000)]"
            + ", [(0:6000000), (0:9000000)]"
            + ", [(0:9000000), (1:2000000)]"
            + ", [(1:2000000), (1:5000000)]"
            + ", [(1:5000000), (1:8000000)]"
            + ", [(1:8000000), (1:11000000)]"
            + ", [(1:11000000), (1:14000000)]"
            + ", [(1:14000000), (1:17000000)]"
            + ", [(1:17000000), (1:20000000)]"
            + ", [(1:20000000), (2:3000000)]"
            + ", [(2:3000000), (2:6000000)]"
            + ", [(2:6000000), (2:9000000)]"
            + ", [(2:9000000), (2:12000000)]"
            + ", [(2:12000000), (2:15000000)]"
            + ", [(2:15000000), (2:18000000)]"
            + ", [(2:18000000), (2:21000000)]"
            + ", [(2:21000000), (2:24000000)]"
            + ", [(2:24000000), (2:27000000)]"
            + ", [(2:27000000), (2:30000000)]]";
    final HashingRegion[] actual;
    actual = HashingRegion.splitWorkload(r, Sex.EITHER, 0, 3, 20, 10000, 30);
    checkRangeArray(actual);
    assertEquals(expected, Arrays.toString(actual));
  }

  static final String REFERENCE_GENOME = ""
          + "#comment" + LS
          + "version" + TAB + "0" + LS
          + "either" + TAB + "def" + TAB + "diploid" + TAB + "linear" + LS
          + "female" + TAB + "seq" + TAB + "seq2" + TAB + "diploid" + TAB + "linear" + LS
          + "male" + TAB + "seq" + TAB + "seq2" + TAB + "none" + TAB + "linear" + LS
          + "";
  public void testNewSplitWorkloadReference() throws IOException {
    final Reader reader = new StringReader(REFERENCE_GENOME);
    final MockSequencesReader r = new MockSequencesReader(SequenceType.DNA, 3);
    final ReferenceGenome rg = new ReferenceGenome(r, reader, Sex.MALE);
    r.setLengths(new int[] {10000000, 20000000, 30000000});
    final String expected = ""
            + "[[(0:0), (0:1500000)]"
            + ", [(0:1500000), (0:3000000)]"
            + ", [(0:3000000), (0:4500000)]"
            + ", [(0:4500000), (0:6000000)]"
            + ", [(0:6000000), (0:7500000)]"
            + ", [(0:7500000), (0:9000000)]"
            + ", [(0:9000000), (1:500000)]"
            + ", [(1:500000), (1:2000000)]"
            + ", [(1:2000000), (1:3500000)]"
            + ", [(1:3500000), (1:5000000)]"
            + ", [(1:5000000), (1:6500000)]"
            + ", [(1:6500000), (1:8000000)]"
            + ", [(1:8000000), (1:9500000)]"
            + ", [(1:9500000), (1:11000000)]"
            + ", [(1:11000000), (1:12500000)]"
            + ", [(1:12500000), (1:14000000)]"
            + ", [(1:14000000), (1:15500000)]"
            + ", [(1:15500000), (1:17000000)]"
            + ", [(1:17000000), (1:18500000)]"
            + ", [(1:18500000), (1:20000000)]]";
    final HashingRegion[] actual;
    actual = HashingRegion.splitWorkload(r, rg, 0, 3, 20, 10000, 30);
    assertEquals(0, actual[0].getStartPaddedPosition());
    assertEquals(1500030, actual[0].getEndPaddedPosition());

    assertEquals(16999970, actual[18].getStartPaddedPosition());
    assertEquals(18500030, actual[18].getEndPaddedPosition());
    checkRangeArray(actual);
    assertEquals(expected, Arrays.toString(actual));
  }

  public void testNewSplitWorkloadSingle() throws IOException {
    final MockSequencesReader r = new MockSequencesReader(SequenceType.DNA, 1);
    r.setLengths(new int[] {10000});
    final String expected = "[[(0:0), (0:10000)]"
        + "]";
    final HashingRegion[] actual;
    actual = HashingRegion.splitWorkload(r, Sex.EITHER, 0, 1, 20, 20000, 30);
    assertEquals(expected, Arrays.toString(actual));
    assertEquals(0, actual[0].getStartPaddedPosition());
    assertEquals(10000, actual[0].getEndPaddedPosition());
  }

  public void testNewSplitWorkloadSingleRange() throws IOException {
    final MockSequencesReader r = new MockSequencesReader(SequenceType.DNA, 3);
    r.setLengths(new int[] {500000, 10000, 600000});
    final String expected = "[[(1:0), (1:10000)]"
        + "]";
    final HashingRegion[] actual;
    actual = HashingRegion.splitWorkload(r, Sex.EITHER, 1, 2, 20, 20000, 30);
    assertEquals(expected, Arrays.toString(actual));
    assertEquals(0, actual[0].getStartPaddedPosition());
    checkRangeArray(actual);
    assertEquals(10000, actual[0].getEndPaddedPosition());
  }


  public void testPaddedInRange() {
    final HashingRegion r = new HashingRegion(4, 5, 5, 6, 3, 8);

    assertEquals(0, r.isInPaddedRange(4, 3));
    assertEquals(-1, r.isInPaddedRange(4, 2));

    assertEquals(0, r.isInPaddedRange(5, 7));
    assertEquals(1, r.isInPaddedRange(5, 8));

    assertEquals(0, r.isInPaddedRange(4, 7));
    assertEquals(0, r.isInPaddedRange(5, 2));

    assertEquals(-1, r.isInPaddedRange(3, 5));
    assertEquals(1, r.isInPaddedRange(6, 5));

    assertEquals(0, HashingRegion.NONE.isInPaddedRange(Long.MAX_VALUE, Long.MAX_VALUE));
    assertEquals(0, HashingRegion.NONE.isInPaddedRange(0, 0));

    final HashingRegion r2 = new HashingRegion(4, 5);
    assertEquals(-1, r2.isInPaddedRange(3, 5));
    assertEquals(1, r2.isInPaddedRange(6, 5));
    assertEquals(0, r2.isInPaddedRange(4, Long.MIN_VALUE));
    assertEquals(1, r2.isInPaddedRange(5, 0));

    final HashingRegion r3 = new HashingRegion(4, 1, 5, 6, 1, 8);
    assertEquals(-1, r3.isInPaddedRange(4, -1));

    final HashingRegion r4 = new HashingRegion(4, 0, 5, 6, 0, 8);
    assertEquals(0, r4.isInPaddedRange(4, -1));
  }
  private static final String PSUEDO_REF = ""
      + "version 1" + LS
      + "either\tdef\thaploid\tlinear" + LS
      + "male\tdup\tseq4:200-400\tseq0:900-1100" + LS
      + "male\tdup\tseq3:200-400\tseq4:550-750" + LS
      ;
  public void testExclusionsPartial() throws IOException {
    final HashingRegion[] regions = {
        new HashingRegion(0, 0, 0, 1000, 0, 1010)
        , new HashingRegion(0, 1001, 1, 500, 990, 510)
        , new HashingRegion(1, 500, 4, 222, 491, 230)
        , new HashingRegion(4, 222, 4, 600, 213, 610)
        , new HashingRegion(4, 600, 4, 700, 591, 710)
        , new HashingRegion(4, 700, 4, 1999, 690, 1999)
    };
    final StringReader reader = new StringReader(PSUEDO_REF);
    final MockSequencesReader genome = new MockSequencesReader(SequenceType.DNA, 5, 10000);
    final ReferenceGenome rg = new ReferenceGenome(genome, reader, Sex.MALE);
    final Map<String, Long> nameMap = new HashMap<>();
    for (long i = 0; i < 5; ++i) {
      nameMap.put("seq" + i, i);
    }
    final List<HashingRegion> excluded = HashingRegion.excludeDuplicateRegions(rg, regions, nameMap);
    // Padded positions aren't checked by equals
    assertEquals(Arrays.asList(
        new HashingRegion(0, 0, 0, 899, 0, 0)
        , new HashingRegion(0, 1100, 1, 500, 0, 0)
        , new HashingRegion(1, 500, 4, 222, 0, 0)
        , new HashingRegion(4, 222, 4, 549, 0, 0)
        , new HashingRegion(4, 750, 4, 1999, 0, 0))
        , excluded);
    assertEquals(899, excluded.get(0).getEndPaddedPosition());
    assertEquals(1100, excluded.get(1).getStartPaddedPosition());
  }

  private static final String PSEUDO_REF_2 = "version 1\n"
          + "either\tdef\tdiploid\tlinear\n"
          + "male\tseq\tseq0\thaploid\tlinear\tseq1\n"
          + "male\tseq\tseq1\thaploid\tlinear\tseq0\n"
          + "male\tdup\tseq0:60000-80000\tseq1:60000-80000\n";

  public void testExclusionCases() throws IOException {
    final StringReader reader = new StringReader(PSEUDO_REF_2);
    final MockSequencesReader genome = new MockSequencesReader(SequenceType.DNA, 3, 300000);
    final ReferenceGenome rg = new ReferenceGenome(genome, reader, Sex.MALE);
    final Map<String, Long> sequenceNameMap = ReaderUtils.getSequenceNameMap(genome);
    //front overlap
    checkExclusionAndPadding(rg, sequenceNameMap, new HashingRegion[]{new HashingRegion(1, 5000, 1, 59999, 4500, 59999)}, new HashingRegion[]{new HashingRegion(1, 5000, 1, 65000, 4500, 65500)});
    //rear overlap
    checkExclusionAndPadding(rg, sequenceNameMap, new HashingRegion[] {new HashingRegion(1, 80000, 1, 90000, 80000, 90500)}, new HashingRegion[] {new HashingRegion(1, 70000, 1, 90000, 69500, 90500)});
    //region entirely within PAR region
    checkExclusionAndPadding(rg, sequenceNameMap, new HashingRegion[] {}, new HashingRegion[] {new HashingRegion(1, 70000, 1, 75000, 69500, 75500)});
    //region entirely within PAR region
    checkExclusionAndPadding(rg, sequenceNameMap, new HashingRegion[] {}, new HashingRegion[] {new HashingRegion(1, 60000, 1, 80000, 59500, 80500)});
    //region entirely within PAR region
    checkExclusionAndPadding(rg, sequenceNameMap, new HashingRegion[] {}, new HashingRegion[] {new HashingRegion(1, 59999, 1, 80000, 59500, 80500)});
    //region entirely consumes PAR region
    checkExclusionAndPadding(rg, sequenceNameMap, new HashingRegion[] {new HashingRegion(1, 55000, 1, 59999, 54500, 59999), new HashingRegion(1, 80000, 1, 85000, 80000, 85500)}, new HashingRegion[] {new HashingRegion(1, 55000, 1, 85000, 54500, 85500)});
    //region on edges
    checkExclusionAndPadding(rg, sequenceNameMap, new HashingRegion[] {
              new HashingRegion(1, 55000, 1, 59999, 54500, 59999),
              new HashingRegion(1, 80000, 1, 85000, 80000, 85500)
            }, new HashingRegion[] {
                    new HashingRegion(1, 55000, 1, 60000, 54500, 60500),
                    new HashingRegion(1, 60000, 1, 80000, 59500, 80500),
                    new HashingRegion(1, 80000, 1, 85000, 79500, 85500)});

    //front overlap padding
    checkExclusionAndPadding(rg, sequenceNameMap, new HashingRegion[]{new HashingRegion(1, 5000, 1, 59999, 4500, 59999)}, new HashingRegion[]{new HashingRegion(1, 5000, 1, 65000, 4500, 70000)});
    //rear overlap padding
    checkExclusionAndPadding(rg, sequenceNameMap, new HashingRegion[]{new HashingRegion(1, 80000, 1, 90000, 80000, 95000)}, new HashingRegion[]{new HashingRegion(1, 70000, 1, 90000, 65000, 95000)});

    checkExclusionAndPadding(rg, sequenceNameMap, new HashingRegion[] {new HashingRegion(0, 50000, 1, 59999, 49500, 59999)}, new HashingRegion[] {new HashingRegion(0, 50000, 1, 65000, 49500, 65500)});
    checkExclusionAndPadding(rg, sequenceNameMap, new HashingRegion[] {new HashingRegion(1, 80000, 2, 50000, 80000, 50500)}, new HashingRegion[] {new HashingRegion(1, 70000, 2, 50000, 69500, 50500)});
    checkExclusionAndPadding(rg, sequenceNameMap, new HashingRegion[] {
            new HashingRegion(0, 50000, 1, 59999, 49500, 59999),
            new HashingRegion(1, 80000, 1, 90000, 80000, 90500)
    }, new HashingRegion[] {new HashingRegion(0, 50000, 1, 90000, 49500, 90500)});
    checkExclusionAndPadding(rg, sequenceNameMap, new HashingRegion[] {
            new HashingRegion(1, 50000, 1, 59999, 49500, 59999),
            new HashingRegion(1, 80000, 2, 50000, 80000, 50500)
    }, new HashingRegion[] {new HashingRegion(1, 50000, 2, 50000, 49500, 50500)});
    checkExclusionAndPadding(rg, sequenceNameMap, new HashingRegion[] {
            new HashingRegion(0, 50000, 1, 59999, 49500, 59999),
            new HashingRegion(1, 80000, 2, 50000, 80000, 50500)
    }, new HashingRegion[] {new HashingRegion(0, 50000, 2, 50000, 49500, 50500)});

    checkExclusionAndPadding(rg, sequenceNameMap, new HashingRegion[] {new HashingRegion(0, 50000, 1, 59999, 49500, 59999)}, new HashingRegion[] {new HashingRegion(0, 50000, 1, 70000, 49500, 70500)});
    checkExclusionAndPadding(rg, sequenceNameMap, new HashingRegion[] {new HashingRegion(1, 80000, 2, 50000, 80000, 50500)}, new HashingRegion[] {new HashingRegion(1, 70000, 2, 50000, 69500, 50500)});

  }

  void checkExclusionAndPadding(ReferenceGenome rg, Map<String, Long> sequenceNameMap, HashingRegion[] expRegion, HashingRegion[] inputRegion) {
    final List<HashingRegion> excluded = HashingRegion.excludeDuplicateRegions(rg, inputRegion, sequenceNameMap);
    assertEquals(Arrays.asList(expRegion), excluded);
    for (int i = 0; i < excluded.size(); ++i) {
      assertEquals(expRegion[i].getStartPaddedPosition(), excluded.get(i).getStartPaddedPosition());
      assertEquals(expRegion[i].getEndPaddedPosition(), excluded.get(i).getEndPaddedPosition());
    }
  }
}


