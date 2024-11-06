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

package com.rtg.variant.bayes.complex;

import com.rtg.reference.Ploidy;
import com.rtg.variant.Variant;
import com.rtg.variant.VariantLocus;
import com.rtg.variant.VariantSample;

import junit.framework.TestCase;

/**
 */
public class IonTorrentCallFilterTest extends TestCase {

  public void testIonFilter() {
    assertFalse(IonTorrentCallFilter.ionFilter(0, 0, 0));
    assertFalse(IonTorrentCallFilter.ionFilter(5, 0, 0));
    assertFalse(IonTorrentCallFilter.ionFilter(0, 5, 0));
    assertFalse(IonTorrentCallFilter.ionFilter(0, 0, 5));

    assertTrue(IonTorrentCallFilter.ionFilter(1, 2, 1));
    assertTrue(IonTorrentCallFilter.ionFilter(2, 4, 3));
  }

  public void testIonTorrentFilter1() {
    final VariantLocus locus = new VariantLocus("foo", 2, 4, "CC", 'N'); //CCC->CC
    final Variant v = new Variant(locus, new VariantSample(Ploidy.DIPLOID, "C:C", false, null, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0));
    assertTrue(IonTorrentCallFilter.ionTorrentFilter(new byte[] {1, 2, 2, 2, 3}, v, true));
  }

  public void testIonTorrentFilter2() {
    final VariantLocus locus = new VariantLocus("foo", 1, 2, "C", 'N'); //CC->C
    final Variant v = new Variant(locus, new VariantSample(Ploidy.DIPLOID, ":", false, null, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0));
    assertTrue(IonTorrentCallFilter.ionTorrentFilter(new byte[] {1, 2, 2, 2, 3}, v, true));
  }

  public void testIonTorrentFilterI() {
    final VariantLocus locus = new VariantLocus("foo", 1, 2, "G", 'N');  //GG->G
    final Variant v = new Variant(locus, new VariantSample(Ploidy.DIPLOID, ":", false, null, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0));
    assertTrue(IonTorrentCallFilter.ionTorrentFilter(new byte[] {4, 3, 3, 4}, v, true));
  }
  public void testIonTorrentFilterI2() {
    final VariantLocus locus = new VariantLocus("foo", 1, 1, "", 'N'); // GG->i
    final Variant v = new Variant(locus, new VariantSample(Ploidy.DIPLOID, ":", false, null, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0));
    assertFalse(IonTorrentCallFilter.ionTorrentFilter(new byte[] {4, 3, 3, 4}, v, true));
  }

  public void testDelete() {
    final VariantLocus locus = new VariantLocus("foo", 4, 5, "C", 'C'); // CCC->CC
    final Variant v = new Variant(locus, new VariantSample(Ploidy.DIPLOID, ":C", false, null, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0));
    assertTrue(IonTorrentCallFilter.ionTorrentFilter(new byte[] {4, 3, 2, 2, 2, 1, 4}, v, true));
  }
  public void testDelete2() {
    final VariantLocus locus = new VariantLocus("foo", 4, 5, "C", 'C'); // CCC->CC
    final Variant v = new Variant(locus, new VariantSample(Ploidy.DIPLOID, "C:", false, null, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0));
    assertTrue(IonTorrentCallFilter.ionTorrentFilter(new byte[] {4, 3, 2, 2, 2, 1, 4}, v, true));
  }
}
