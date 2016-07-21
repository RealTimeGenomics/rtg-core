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
    final IndelMatcher m = new IndelMatcher(new ReferenceBasedBuffer<>(1 - 0, IndelDetectorFactory.SINGLETON, template, 0));
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
    final IndelMatcher m = new IndelMatcher(new ReferenceBasedBuffer<>(1 - 0, IndelDetectorFactory.SINGLETON, template, 0));
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
      for (int i = 0; i < count; i++) {
        m.match(pos, evid);
      }
      pos++;
    }
  }
}

