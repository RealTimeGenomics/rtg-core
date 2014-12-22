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
package com.rtg.calibrate;

import java.util.Arrays;

import com.rtg.util.TestUtils;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMSequenceRecord;

import junit.framework.TestCase;

/**
 */
public class CovariateEnumTest extends TestCase {

  public void test() {
    TestUtils.testEnum(CovariateEnum.class, "[READGROUP, READPOSITION, BASEQUALITY, SEQUENCE]");
  }

  public void testDefaults() {
    assertEquals(3, CovariateEnum.DEFAULT_COVARIATES.size());
    assertEquals(CovariateEnum.READGROUP, CovariateEnum.DEFAULT_COVARIATES.get(0));
    assertEquals(CovariateEnum.SEQUENCE, CovariateEnum.DEFAULT_COVARIATES.get(1));
    assertEquals(CovariateEnum.BASEQUALITY, CovariateEnum.DEFAULT_COVARIATES.get(2));
  }

  public void testNormal() {
    final Covariate[] covs = CovariateEnum.getCovariates(Arrays.asList(CovariateEnum.SEQUENCE, CovariateEnum.READGROUP, CovariateEnum.READPOSITION, CovariateEnum.READPOSITION, CovariateEnum.BASEQUALITY), null);
    assertEquals(4, covs.length);
    for (int i = 0; i < covs.length; i++) {
      assertEquals(CovariateEnum.values()[i], covs[i].getType());
    }
    assertEquals(0, covs[CovariateEnum.READPOSITION.ordinal()].size());
  }

  public void testWithHeader() {
    final SAMFileHeader header = new SAMFileHeader();
    header.addReadGroup(new SAMReadGroupRecord("readGroup1"));
    header.getSequenceDictionary().addSequence(new SAMSequenceRecord("sequence1", 1));
    Covariate[] covs = CovariateEnum.getCovariates(CovariateEnum.DEFAULT_COVARIATES, header);
    assertEquals(3, covs.length);
    assertEquals(CovariateEnum.READGROUP, covs[0].getType());
    assertTrue(covs[0] instanceof CovariateSingleReadGroup);
    assertEquals(CovariateEnum.BASEQUALITY, covs[1].getType());
    assertEquals(CovariateEnum.SEQUENCE, covs[2].getType());
    assertTrue(covs[2] instanceof CovariateSequenceFixed);
    header.addReadGroup(new SAMReadGroupRecord("readGroup2"));
    covs = CovariateEnum.getCovariates(Arrays.asList(CovariateEnum.READGROUP), header);
    assertEquals(1, covs.length);
    assertEquals(CovariateEnum.READGROUP, covs[0].getType());
    assertTrue(covs[0] instanceof CovariateReadGroup);
  }
}
