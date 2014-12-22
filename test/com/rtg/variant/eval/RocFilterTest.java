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

package com.rtg.variant.eval;

import com.rtg.util.TestUtils;
import com.rtg.vcf.VcfReader;
import com.rtg.vcf.VcfRecord;

import junit.framework.TestCase;

/**
 */
public class RocFilterTest extends TestCase {

  private static final VcfRecord PASS_RECORD = VcfReader.vcfLineToRecord("chr1 250 . A C 20.0 PASS . GT 1/1".replaceAll(" ", "\t"));
  private static final VcfRecord FAIL_RECORD = VcfReader.vcfLineToRecord("chr1 250 . A C 20.0 RC . GT 1/1".replaceAll(" ", "\t"));

  private static final VcfRecord HOMOZYGOUS_RECORD = VcfReader.vcfLineToRecord("chr1 250 . A C 20.0 PASS . GT 1/1".replaceAll(" ", "\t"));
  private static final VcfRecord HETEROZYGOUS_RECORD = VcfReader.vcfLineToRecord("chr1 250 . A C 20.0 PASS . GT 0/1".replaceAll(" ", "\t"));
  private static final VcfRecord IDENTITY_RECORD = VcfReader.vcfLineToRecord("chr1 250 . A . 20.0 PASS . GT 0/0".replaceAll(" ", "\t"));

  private static final VcfRecord COMPLEX_HOMOZYGOUS_RECORD = VcfReader.vcfLineToRecord("chr1 250 . A C 20.0 PASS XRX GT 1/1".replaceAll(" ", "\t"));
  private static final VcfRecord COMPLEX_HETEROZYGOUS_RECORD = VcfReader.vcfLineToRecord("chr1 250 . A C 20.0 PASS XRX GT 0/1".replaceAll(" ", "\t"));
  private static final VcfRecord COMPLEX_IDENTITY_RECORD = VcfReader.vcfLineToRecord("chr1 250 . A . 20.0 PASS XRX GT 0/0".replaceAll(" ", "\t"));

  public void testEnum() {
    TestUtils.testEnum(RocFilter.class, "[ALL, HOMOZYGOUS, HETEROZYGOUS, COMPLEX, SIMPLE, HOMOZYGOUS_COMPLEX, HOMOZYGOUS_SIMPLE, HETEROZYGOUS_COMPLEX, HETEROZYGOUS_SIMPLE]");
  }

  public void testAll() {
    final RocFilter f = RocFilter.ALL;
    assertTrue(f.accept(PASS_RECORD, 0));
    assertTrue(f.accept(FAIL_RECORD, 0));
    assertTrue(f.accept(HOMOZYGOUS_RECORD, 0));
    assertTrue(f.accept(HETEROZYGOUS_RECORD, 0));
    assertTrue(f.accept(IDENTITY_RECORD, 0));

    assertTrue(f.accept(COMPLEX_HETEROZYGOUS_RECORD, 0));
    assertTrue(f.accept(COMPLEX_HOMOZYGOUS_RECORD, 0));
    assertTrue(f.accept(COMPLEX_IDENTITY_RECORD, 0));
  }

  public void testHomozygous() {
    final RocFilter f = RocFilter.HOMOZYGOUS;
    assertTrue(f.accept(PASS_RECORD, 0));
    assertTrue(f.accept(FAIL_RECORD, 0));
    assertTrue(f.accept(HOMOZYGOUS_RECORD, 0));
    assertTrue(f.accept(COMPLEX_HOMOZYGOUS_RECORD, 0));

    assertFalse(f.accept(HETEROZYGOUS_RECORD, 0));
    assertFalse(f.accept(IDENTITY_RECORD, 0));
    assertFalse(f.accept(COMPLEX_HETEROZYGOUS_RECORD, 0));
    assertFalse(f.accept(COMPLEX_IDENTITY_RECORD, 0));
  }

  public void testHeterozygous() {
    final RocFilter f = RocFilter.HETEROZYGOUS;
    assertTrue(f.accept(HETEROZYGOUS_RECORD, 0));
    assertTrue(f.accept(COMPLEX_HETEROZYGOUS_RECORD, 0));

    assertFalse(f.accept(PASS_RECORD, 0));
    assertFalse(f.accept(FAIL_RECORD, 0));
    assertFalse(f.accept(HOMOZYGOUS_RECORD, 0));
    assertFalse(f.accept(COMPLEX_HOMOZYGOUS_RECORD, 0));
    assertFalse(f.accept(IDENTITY_RECORD, 0));
    assertFalse(f.accept(COMPLEX_IDENTITY_RECORD, 0));
  }

  public void testComplex() {
    final RocFilter f = RocFilter.COMPLEX;
    assertTrue(f.accept(COMPLEX_HETEROZYGOUS_RECORD, 0));
    assertTrue(f.accept(COMPLEX_HOMOZYGOUS_RECORD, 0));
    assertTrue(f.accept(COMPLEX_IDENTITY_RECORD, 0));

    assertFalse(f.accept(HETEROZYGOUS_RECORD, 0));
    assertFalse(f.accept(PASS_RECORD, 0));
    assertFalse(f.accept(FAIL_RECORD, 0));
    assertFalse(f.accept(HOMOZYGOUS_RECORD, 0));
    assertFalse(f.accept(IDENTITY_RECORD, 0));
  }


  public void testSimple() {
    final RocFilter f = RocFilter.SIMPLE;
    assertTrue(f.accept(HETEROZYGOUS_RECORD, 0));
    assertTrue(f.accept(FAIL_RECORD, 0));
    assertTrue(f.accept(PASS_RECORD, 0));
    assertTrue(f.accept(HOMOZYGOUS_RECORD, 0));
    assertTrue(f.accept(IDENTITY_RECORD, 0));

    assertFalse(f.accept(COMPLEX_HETEROZYGOUS_RECORD, 0));
    assertFalse(f.accept(COMPLEX_HOMOZYGOUS_RECORD, 0));
    assertFalse(f.accept(COMPLEX_IDENTITY_RECORD, 0));
  }

  public void testHomozygousSimple() {
    final RocFilter f = RocFilter.HOMOZYGOUS_SIMPLE;
    assertTrue(f.accept(FAIL_RECORD, 0));
    assertTrue(f.accept(PASS_RECORD, 0));
    assertTrue(f.accept(HOMOZYGOUS_RECORD, 0));

    assertFalse(f.accept(HETEROZYGOUS_RECORD, 0));
    assertFalse(f.accept(COMPLEX_HETEROZYGOUS_RECORD, 0));
    assertFalse(f.accept(COMPLEX_HOMOZYGOUS_RECORD, 0));
    assertFalse(f.accept(COMPLEX_IDENTITY_RECORD, 0));
    assertFalse(f.accept(IDENTITY_RECORD, 0));
  }

  public void testHomozygousComplex() {
    final RocFilter f = RocFilter.HOMOZYGOUS_COMPLEX;
    assertTrue(f.accept(COMPLEX_HOMOZYGOUS_RECORD, 0));

    assertFalse(f.accept(FAIL_RECORD, 0));
    assertFalse(f.accept(PASS_RECORD, 0));
    assertFalse(f.accept(HOMOZYGOUS_RECORD, 0));
    assertFalse(f.accept(HETEROZYGOUS_RECORD, 0));
    assertFalse(f.accept(COMPLEX_HETEROZYGOUS_RECORD, 0));
    assertFalse(f.accept(COMPLEX_IDENTITY_RECORD, 0));
    assertFalse(f.accept(IDENTITY_RECORD, 0));
  }

  public void testHeterozygousSimple() {
    final RocFilter f = RocFilter.HETEROZYGOUS_SIMPLE;
    assertTrue(f.accept(HETEROZYGOUS_RECORD, 0));

    assertFalse(f.accept(FAIL_RECORD, 0));
    assertFalse(f.accept(PASS_RECORD, 0));
    assertFalse(f.accept(HOMOZYGOUS_RECORD, 0));
    assertFalse(f.accept(COMPLEX_HETEROZYGOUS_RECORD, 0));
    assertFalse(f.accept(COMPLEX_HOMOZYGOUS_RECORD, 0));
    assertFalse(f.accept(COMPLEX_IDENTITY_RECORD, 0));
    assertFalse(f.accept(IDENTITY_RECORD, 0));
  }

  public void testHeterozygousComplex() {
    final RocFilter f = RocFilter.HETEROZYGOUS_COMPLEX;
    assertTrue(f.accept(COMPLEX_HETEROZYGOUS_RECORD, 0));

    assertFalse(f.accept(FAIL_RECORD, 0));
    assertFalse(f.accept(PASS_RECORD, 0));
    assertFalse(f.accept(HOMOZYGOUS_RECORD, 0));
    assertFalse(f.accept(HETEROZYGOUS_RECORD, 0));
    assertFalse(f.accept(COMPLEX_HOMOZYGOUS_RECORD, 0));
    assertFalse(f.accept(COMPLEX_IDENTITY_RECORD, 0));
    assertFalse(f.accept(IDENTITY_RECORD, 0));
  }

}
