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

package com.rtg.variant.bayes.multisample.cancer;

import static com.rtg.variant.bayes.multisample.cancer.SomaticFilterTest.getSomaticVcfRecord;

import com.rtg.vcf.VcfRecord;

import junit.framework.TestCase;

/**
 */
public class SomaticRecordUtilsTest extends TestCase {

  public void testLohRefRef() {
    final VcfRecord record = getSomaticVcfRecord("0/0", "0/0");
    assertEquals(0.0, SomaticRecordUtils.lossOfHeterozygosity(record));
  }

  public void testLohRefAlt() {
    final VcfRecord record = getSomaticVcfRecord("0/0", "0/1");
    assertEquals(-1.0, SomaticRecordUtils.lossOfHeterozygosity(record));
  }

  public void testLohAltRef() {
    final VcfRecord record = getSomaticVcfRecord("0/1", "0/0");
    assertEquals(1.0, SomaticRecordUtils.lossOfHeterozygosity(record));
  }

  public void testLohAltAlt() {
    final VcfRecord record = getSomaticVcfRecord("0/1", "0/1");
    assertEquals(0.0, SomaticRecordUtils.lossOfHeterozygosity(record));
  }

  public void testLohAltAltSwapped() {
    final VcfRecord record = getSomaticVcfRecord("0/1", "1/0");
    assertEquals(0.0, SomaticRecordUtils.lossOfHeterozygosity(record));
  }

  public void testLohHaploidAltAlt() {
    final VcfRecord record = getSomaticVcfRecord("1", "1");
    assertEquals(0.0, SomaticRecordUtils.lossOfHeterozygosity(record));
  }

  public void testLohHaploidAltRef() {
    final VcfRecord record = getSomaticVcfRecord("1", "0");
    assertEquals(-1.0, SomaticRecordUtils.lossOfHeterozygosity(record));
  }

  public void testIsGainOfRef() {
    final VcfRecord record = getSomaticVcfRecord("0/1", "0/0");
    assertTrue(SomaticRecordUtils.isGainOfReference(record));
  }

  public void testIsNotGainOfRef() {
    final VcfRecord record = getSomaticVcfRecord("0/1", "1/1");
    assertFalse(SomaticRecordUtils.isGainOfReference(record));
  }

}