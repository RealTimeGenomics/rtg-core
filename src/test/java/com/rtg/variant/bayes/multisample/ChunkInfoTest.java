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

package com.rtg.variant.bayes.multisample;

import com.rtg.sam.SamFilterParams;
import com.rtg.variant.VariantParams;

import junit.framework.TestCase;

/**
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
