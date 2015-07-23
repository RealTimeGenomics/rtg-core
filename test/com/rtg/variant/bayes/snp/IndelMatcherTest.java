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
    final IndelMatcher m = new IndelMatcher(template, 0, 1);
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
}
