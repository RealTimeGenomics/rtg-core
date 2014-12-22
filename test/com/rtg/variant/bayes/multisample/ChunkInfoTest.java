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

package com.rtg.variant.bayes.multisample;

import com.rtg.sam.SamFilterParams;
import com.rtg.variant.VariantParams;

import junit.framework.TestCase;

/**
 *         Date: 7/03/12
 *         Time: 4:42 PM
 */
public class ChunkInfoTest extends TestCase {
  public void testConstructor() {
    final VariantParams vp = VariantParams.builder().create();
    final ChunkInfo chunks = new ChunkInfo(100, "foo", 10, vp.filterParams().restrictionStart(), vp.filterParams().restrictionEnd(), 5, 4);
    assertEquals(0, chunks.start());
    assertEquals(100, chunks.end());
    assertEquals(100, chunks.length());
    assertEquals(10, chunks.chunkSize());
    assertEquals(10, chunks.numberChunks());
    assertEquals(50, chunks.percent(50));
    assertEquals(144, chunks.bufferSize());
  }
  public void testRestriction() {
    final VariantParams vp = VariantParams.builder().filterParams(SamFilterParams.builder().restriction("foo:20-80").create()).create();
    final ChunkInfo chunks = new ChunkInfo(100, "foo", 10, vp.filterParams().restrictionStart(), vp.filterParams().restrictionEnd(), 5, 4);
    assertEquals(19, chunks.start());
    assertEquals(80, chunks.end());
    assertEquals(61, chunks.length());
    assertEquals(10, chunks.chunkSize());
    assertEquals(7, chunks.numberChunks());
    assertEquals(50, chunks.percent(50));
    assertEquals(144, chunks.bufferSize());
  }
  public void testRestrictionSimple() {
    final VariantParams vp = VariantParams.builder().filterParams(SamFilterParams.builder().restriction("foo:20-20").create()).create();
    final ChunkInfo chunks = new ChunkInfo(100, "foo", 10, vp.filterParams().restrictionStart(), vp.filterParams().restrictionEnd(), 5, 4);
    assertEquals(19, chunks.start());
    assertEquals(20, chunks.end());
    assertEquals(1, chunks.length());
  }
}
