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

package com.rtg.variant;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import com.rtg.reference.Ploidy;
import com.rtg.util.TestUtils;
import com.rtg.variant.Variant.VariantFilter;
import com.rtg.variant.bayes.ArrayGenotypeMeasure;
import com.rtg.variant.bayes.GenotypeMeasure;
import com.rtg.variant.bayes.snp.HypothesesSnp;
import com.rtg.variant.util.VariantUtils;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 */
public class VariantTest extends TestCase {

  public void testSetters() {
    final Variant v = new Variant(new VariantLocus("foo", 1, 2, "c", (char) -1));
    assertFalse(v.isInteresting());
    v.setInteresting();
    assertTrue(v.isInteresting());

    assertFalse(v.isComplexScored());
    v.setComplexScored();
    assertTrue(v.isComplexScored());

    assertNull(v.getNonIdentityPosterior());
    v.setNonIdentityPosterior(-1.5);
    assertEquals(-1.5, v.getNonIdentityPosterior());

    assertNull(v.getPossibleCause());
    assertNull(v.getPossibleCauseScore());
    v.setPossibleCause("dd");
    v.setPossibleCauseScore(1.5);
    assertEquals("dd", v.getPossibleCause());
    assertEquals(1.5, v.getPossibleCauseScore());

    assertFalse(v.isComplexEquivalent());
    v.setComplexEquivalent();
    assertTrue(v.isComplexEquivalent());

    assertFalse(v.isIndel());
    v.setIndel(0);
    assertTrue(v.isIndel());

    assertFalse(v.isOverflow());
    v.setOverflow();
    assertTrue(v.isOverflow());
  }

  public void testFilters() {
    final Variant v = new Variant(new VariantLocus("foo", 1, 2, "c", (char) -1));
    assertFalse(v.isFiltered());
    assertFalse(v.isFiltered(VariantFilter.AMBIGUITY));
    v.addFilter(VariantFilter.AMBIGUITY);
    assertTrue(v.isFiltered());
    assertTrue(v.isFiltered(VariantFilter.AMBIGUITY));

    final Variant v2 = new Variant(new VariantLocus("foo", 1, 2, "c", (char) -1), VariantFilter.COVERAGE.mask());
    assertTrue(v2.isFiltered());
    assertTrue(v2.isFiltered(VariantFilter.COVERAGE));
    assertFalse(v2.isFiltered(VariantFilter.AMBIGUITY));
  }

  public void testOrder() {
    final Variant v = new Variant(new VariantLocus("foo", 1, 2), new VariantSample(Ploidy.DIPLOID));
    final Variant v2 = new Variant(new VariantLocus("foo", 1, 2), VariantFilter.COVERAGE.mask(), new VariantSample(Ploidy.DIPLOID));
    v2.setIndel(0);
    TestUtils.testOrder(new Variant[] {v2, v}, false);
  }

  public void testCopy() {
    final Variant v = new Variant(new VariantLocus("foo", 1, 2), new VariantSample(Ploidy.DIPLOID));
    v.setComplexEquivalent();
    v.setComplexScored();
    v.setIndel(0);
    v.setInteresting();
    v.setNonIdentityPosterior(12.5);
    v.setPossibleCause("CA");
    v.setPossibleCauseScore(44.5);
    v.setNormalCancerScore(45.6);
    v.setDiseasePresenceScore(46.7);
    v.addFilter(VariantFilter.COVERAGE);

    final Variant v2 = new Variant(new VariantLocus("foo", 1, 2), new VariantSample(Ploidy.DIPLOID));
    assertFalse(v2.isComplexEquivalent());
    assertFalse(v2.isComplexScored());
    assertFalse(v2.isIndel());
    assertFalse(v2.isFiltered());
    assertNull(v2.getNonIdentityPosterior());
    assertNull(v2.getPossibleCause());
    assertNull(v2.getPossibleCauseScore());

    Variant.copy(v, v2);

    assertTrue(v2.isComplexEquivalent());
    assertTrue(v2.isComplexScored());
    assertTrue(v2.isIndel());
    assertTrue(v2.isFiltered());
    assertEquals(12.5, v2.getNonIdentityPosterior());
    assertEquals("CA", v2.getPossibleCause());
    assertEquals(44.5, v2.getPossibleCauseScore());
    assertEquals(45.6, v2.getNormalCancerScore());
    assertEquals(46.7, v2.getDiseasePresenceScore());
    assertEquals(VariantFilter.COVERAGE.mask(), v2.getFilters());
  }

  public void testIsSnp() {
    Variant v = new Variant(new VariantLocus("foo", 1, 3), new VariantSample(Ploidy.DIPLOID, "a", true, null, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0));
    assertFalse(v.isSnp()); //null locus name

    v = new Variant(new VariantLocus("", 1, 1, "", (char) -1), new VariantSample(Ploidy.DIPLOID, "a", true, null, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0));
    assertFalse(v.isSnp()); //0 length locus name

    v = new Variant(new VariantLocus("", 1, 1, "", (char) -1), new VariantSample(Ploidy.DIPLOID, "a", true, null, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0));
    assertFalse(v.isSnp()); //ref is insert

    v = new Variant(new VariantLocus("", 1, 3, "cc", (char) -1), new VariantSample(Ploidy.DIPLOID, "a", true, null, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0));
    assertFalse(v.isSnp());  //distance between start & end too great

    v = new Variant(new VariantLocus("", 1, 2, "c", (char) -1), new VariantSample(Ploidy.DIPLOID, "ag", true, null, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0));
    assertFalse(v.isSnp());  //bc name is mnp

    v = new Variant(new VariantLocus("", 1, 2, "c", (char) -1), new VariantSample(Ploidy.DIPLOID, "agt", true, null, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0));
    assertFalse(v.isSnp());  //bc name is mnp

    v = new Variant(new VariantLocus("", 1, 2, "c", (char) -1), new VariantSample(Ploidy.DIPLOID));
    v.addFilter(VariantFilter.FAILED_COMPLEX);
    assertFalse(v.isSnp()); //failed complex.

    v = new Variant(new VariantLocus("", 1, 2, "c", (char) -1), new VariantSample(Ploidy.DIPLOID));  //no bestCategory, not filtered by failed complex
    assertTrue(v.isSnp());

    v = new Variant(new VariantLocus("", 1, 2, "c", (char) -1), new VariantSample(Ploidy.DIPLOID, "a", true, null, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0));
    assertTrue(v.isSnp());

    v = new Variant(new VariantLocus("", 1, 2, "c", (char) -1), new VariantSample(Ploidy.DIPLOID, "a" + VariantUtils.COLON + "t", false, null, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0));
    assertTrue(v.isSnp());
  }

  public void testAlleleSet() {
    // Was testing disagreeing hypothesis but that doesn't exist now.
    final GenotypeMeasure aMeasure = new ArrayGenotypeMeasure(SimplePossibility.SINGLETON, new double[] {0.3, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1}, new HypothesesSnp(SimplePossibility.SINGLETON, GenomePriorParams.builder().create(), false, 0));
    final GenotypeMeasure bMeasure = new ArrayGenotypeMeasure(SimplePossibility.SINGLETON, new double[] {0.1, 0.1, 0.1, 0.1, 0.3, 0.1, 0.1, 0.1, 0.1, 0.1}, new HypothesesSnp(SimplePossibility.SINGLETON, GenomePriorParams.builder().create(), false, 0));
    final VariantSample a = new VariantSample(Ploidy.DIPLOID, "A:A", true, aMeasure, VariantSample.DeNovoStatus.NOT_DE_NOVO, null);
    final VariantSample b = new VariantSample(Ploidy.DIPLOID, "A:C", true, bMeasure, VariantSample.DeNovoStatus.NOT_DE_NOVO, null);
    final Set<String> alleles = Variant.alleleSet(new VariantSample[]{a, b}, "A");
    final Set<String> expected = new HashSet<>(Arrays.asList("A", "C"));
    assertEquals(expected, alleles);
  }
  public void testAlleleSetNoRef() {
    final GenotypeMeasure aMeasure = new ArrayGenotypeMeasure(SimplePossibility.SINGLETON, new double[] {0.3, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1}, new HypothesesSnp(SimplePossibility.SINGLETON, GenomePriorParams.builder().create(), false, 3));
    final GenotypeMeasure bMeasure = new ArrayGenotypeMeasure(SimplePossibility.SINGLETON, new double[] {0.3, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1}, new HypothesesSnp(SimplePossibility.SINGLETON, GenomePriorParams.builder().create(), false, 3));
    final VariantSample a = new VariantSample(Ploidy.DIPLOID, "A:A", true, aMeasure, VariantSample.DeNovoStatus.NOT_DE_NOVO, null);
    final VariantSample b = new VariantSample(Ploidy.DIPLOID, "A:A", true, bMeasure, VariantSample.DeNovoStatus.NOT_DE_NOVO, null);
    final Set<String> alleles = Variant.alleleSet(new VariantSample[]{a, b}, "T");
    final Set<String> expected = new HashSet<>(Arrays.asList("A", "T"));
    assertEquals(expected, alleles);
  }
  public void testAlleleSetNRef() {
    final GenotypeMeasure aMeasure = new ArrayGenotypeMeasure(SimplePossibility.SINGLETON, new double[] {0.3, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1}, new HypothesesSnp(SimplePossibility.SINGLETON, GenomePriorParams.builder().create(), false, -1));
    final GenotypeMeasure bMeasure = new ArrayGenotypeMeasure(SimplePossibility.SINGLETON, new double[] {0.3, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1}, new HypothesesSnp(SimplePossibility.SINGLETON, GenomePriorParams.builder().create(), false, -1));
    final VariantSample a = new VariantSample(Ploidy.DIPLOID, "A:A", true, aMeasure, VariantSample.DeNovoStatus.NOT_DE_NOVO, null);
    final VariantSample b = new VariantSample(Ploidy.DIPLOID, "A:A", true, bMeasure, VariantSample.DeNovoStatus.NOT_DE_NOVO, null);
    final Set<String> alleles = Variant.alleleSet(new VariantSample[]{a, b}, "N");
    final Set<String> expected = new HashSet<>(Arrays.asList("A", "N"));
    assertEquals(expected, alleles);
  }
}
