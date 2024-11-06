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
    final EvidenceInterface di = new EvidenceQ(DescriptionSnp.SINGLETON, 0, 0, 0, 0.1, 0.1, true, false, true, false, false);
    mo.increment(di);
  }

  public void testDiploid() {
    final GenomePriorParams params = GenomePriorParams.builder().create();

    final ModelSnpFactory mf = new ModelSnpFactory(params, false, new NoAlleleBalance());
    mf.globalIntegrity();
    final ModelInterface<?> mo = mf.make(0);
    assertEquals(10, mo.size());
    final EvidenceInterface di = new EvidenceQ(DescriptionSnp.SINGLETON, 0, 0, 0, 0.1, 0.1, true, false, true, false, false);
    mo.increment(di);
  }
}
