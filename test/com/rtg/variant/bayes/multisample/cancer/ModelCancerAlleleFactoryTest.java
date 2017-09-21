/*
 * Copyright (c) 2017. Real Time Genomics Limited.
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

import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.EvidenceInterface;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.NoAlleleBalance;
import com.rtg.variant.bayes.snp.DescriptionSnp;
import com.rtg.variant.bayes.snp.EvidenceQ;

import junit.framework.TestCase;

/**
 */
public class ModelCancerAlleleFactoryTest extends TestCase {

  private void checkRef(final int refNt) {
    final GenomePriorParams params = GenomePriorParams.builder().create();
    final ModelCancerAlleleFactory mf = new ModelCancerAlleleFactory(params, false, new NoAlleleBalance());
    mf.globalIntegrity();
    final ModelInterface<Description> mo = mf.make(refNt);
    assertEquals(15, mo.size());
    final EvidenceInterface di = new EvidenceQ(DescriptionSnp.SINGLETON, 0, 0, 0, 0.1, 0.1, true, false, false, false);
    mo.increment(di);
  }

  public void test() {
    checkRef(0);
    checkRef(-1);
  }
}
