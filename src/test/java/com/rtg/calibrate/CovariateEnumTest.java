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
package com.rtg.calibrate;

import java.util.Arrays;
import java.util.Collections;

import com.rtg.util.TestUtils;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMSequenceRecord;

import junit.framework.TestCase;

/**
 */
public class CovariateEnumTest extends TestCase {

  public void test() {
    TestUtils.testEnum(CovariateEnum.class, "[READGROUP, BASEQUALITY, SEQUENCE, ARM, MACHINECYCLE]");
  }

  public void testDefaults() {
    assertEquals(3, CovariateEnum.DEFAULT_COVARIATES.size());
    assertEquals(CovariateEnum.READGROUP, CovariateEnum.DEFAULT_COVARIATES.get(0));
    assertEquals(CovariateEnum.SEQUENCE, CovariateEnum.DEFAULT_COVARIATES.get(1));
    assertEquals(CovariateEnum.BASEQUALITY, CovariateEnum.DEFAULT_COVARIATES.get(2));
  }

  public void testNormal() {
    final Covariate[] covs = CovariateEnum.getCovariates(Arrays.asList(CovariateEnum.SEQUENCE, CovariateEnum.READGROUP, CovariateEnum.ARM, CovariateEnum.ARM, CovariateEnum.BASEQUALITY), null);
    assertEquals(4, covs.length);
    for (int i = 0; i < covs.length; ++i) {
      assertEquals(CovariateEnum.values()[i], covs[i].getType());
    }
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
    covs = CovariateEnum.getCovariates(Collections.singletonList(CovariateEnum.READGROUP), header);
    assertEquals(1, covs.length);
    assertEquals(CovariateEnum.READGROUP, covs[0].getType());
    assertTrue(covs[0] instanceof CovariateReadGroup);
  }
}
