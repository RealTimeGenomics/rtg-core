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
package com.rtg.variant.bayes.snp;

import com.rtg.variant.Variant;
import com.rtg.variant.VariantOutputLevel;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantParamsBuilder;

import junit.framework.TestCase;

/**
 */
public class IndelMatcherTest extends TestCase {

  public void test() {
    final VariantParams params = new VariantParamsBuilder().callLevel(VariantOutputLevel.ALL).create();
    final byte[] template = {1, 2, 3, 4, 0};
    final IndelMatcher m = new IndelMatcher(new ReferenceBasedBuffer<>(1, IndelDetectorFactory.SINGLETON, template, 0));
    assertNotNull(m);
    final double phred = 0.0001;
    m.match(0, new EvidenceIndel(phred, EvidenceIndel.INSERT, 0));
    m.match(0, new EvidenceIndel(phred, EvidenceIndel.INSERT, 0));
    m.match(0, new EvidenceIndel(phred, EvidenceIndel.INSERT, 0));
    m.match(0, new EvidenceIndel(phred, EvidenceIndel.INSERT, 0));
    m.match(0, new EvidenceIndel(phred, EvidenceIndel.INSERT, 0));
    m.match(0, new EvidenceIndel(phred, EvidenceIndel.INSERT, 0));
    m.match(1, null);
    assertNotNull(m.output("foo", 0, params));
  }

  public void testDelete() {
    final Variant foo = getVariant(EvidenceIndel.DELETE);
    assertTrue(foo.isInteresting());
    assertTrue(foo.isIndel());
    assertEquals("foo", foo.getLocus().getSequenceName());
    assertEquals(1, foo.getLocus().getStart());
    assertEquals(2, foo.getLocus().getEnd());
    assertFalse(foo.isSoftClip());
  }

  public void testSoftClipLeft() {
    final Variant foo = getVariant(EvidenceIndel.SOFT_CLIP_LEFT);
    assertTrue(foo.isInteresting());
    assertFalse(foo.isIndel());
    assertEquals("foo", foo.getLocus().getSequenceName());
    assertEquals(1, foo.getLocus().getStart());
    assertEquals(1, foo.getLocus().getEnd());
    assertTrue(foo.isSoftClip());
    assertEquals(Variant.SoftClipSide.LEFT, foo.getSoftClipSide());
  }

  public void testSoftClipRight() {
    final int type = EvidenceIndel.SOFT_CLIP_RIGHT;

    final Variant foo = getVariant(type);
    assertTrue(foo.isInteresting());
    assertFalse(foo.isIndel());
    assertEquals("foo", foo.getLocus().getSequenceName());
    assertEquals(1, foo.getLocus().getStart());
    assertEquals(1, foo.getLocus().getEnd());
    assertTrue(foo.isSoftClip());
    assertEquals(Variant.SoftClipSide.RIGHT, foo.getSoftClipSide());
  }

  private Variant getVariant(int type) {
    final VariantParams params = new VariantParamsBuilder().callLevel(VariantOutputLevel.ALL).create();
    final byte[] template = {1, 2, 3, 4, 4, 0};
    final IndelMatcher m = new IndelMatcher(new ReferenceBasedBuffer<>(1, IndelDetectorFactory.SINGLETON, template, 0));
    assertNotNull(m);
    final double phred = 0.0001;
    m.match(0, null);
    match(m, 1, new int[] {2, 2}, new EvidenceIndel(phred, type, 1));
    m.match(3, null);
    m.match(4, null);
    m.match(5, null);
    m.output("foo", 0, params);
    return m.output("foo", 1, params);
  }

  private static void match(IndelMatcher m, int start, int[] counts, EvidenceIndel evid) {
    int pos = start;
    for (int count : counts) {
      for (int i = 0; i < count; ++i) {
        m.match(pos, evid);
      }
      ++pos;
    }
  }
}

