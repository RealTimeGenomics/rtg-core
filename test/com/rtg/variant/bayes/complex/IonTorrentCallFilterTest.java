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
