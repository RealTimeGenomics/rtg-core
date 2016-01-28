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

import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.bayes.EvidenceInterface;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.NoAlleleBalance;

import junit.framework.TestCase;

/**
 */
public class ModelSnpFactoryTest extends TestCase {

  public void testHaploid() {
    final GenomePriorParams params = GenomePriorParams.builder().create();

    final ModelSnpFactory mf = new ModelSnpFactory(params, true, new NoAlleleBalance());
    mf.globalIntegrity();
    final ModelInterface<?> mo = mf.make(0);
    assertEquals(4, mo.size());
    final EvidenceInterface di = new EvidenceQ(DescriptionSnp.SINGLETON, 0, 0, 0, 0.1, 0.1, true, false, false, false);
    mo.increment(di);
  }

  public void testDiploid() {
    final GenomePriorParams params = GenomePriorParams.builder().create();

    final ModelSnpFactory mf = new ModelSnpFactory(params, false, new NoAlleleBalance());
    mf.globalIntegrity();
    final ModelInterface<?> mo = mf.make(0);
    assertEquals(10, mo.size());
    final EvidenceInterface di = new EvidenceQ(DescriptionSnp.SINGLETON, 0, 0, 0, 0.1, 0.1, true, false, false, false);
    mo.increment(di);
  }
}
