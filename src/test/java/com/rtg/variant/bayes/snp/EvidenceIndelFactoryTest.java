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

import com.rtg.variant.bayes.EvidenceInterface;
import com.rtg.variant.bayes.Model;

import junit.framework.TestCase;

/**
 * Test class
 */
public class EvidenceIndelFactoryTest extends TestCase {
  public void test() {
    final CachedEvidenceFactory fact = EvidenceIndelFactory.SINGLETON;
    final EvidenceInterface ei = fact.evidence(1, 0, 0, 4, 7, 0, 0, false);
    assertEquals(0.0, ei.mapError());
    final EvidenceInterface ei2 = fact.evidence(1, 0, 0, Model.AMBIGUITY_PHRED - 1, 7, 0, 0, false);
    assertEquals(Model.AMBIGUITY_THRESHOLD, ei2.mapError());
  }
}
