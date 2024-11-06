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

package com.rtg.variant.bayes.multisample;


import com.rtg.launcher.MockReaderParams;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reference.SexMemo;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.variant.GenomePriorParamsBuilder;
import com.rtg.variant.MachineErrorChooserInterface;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantParamsBuilder;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.NoAlleleBalance;
import com.rtg.variant.bayes.snp.ModelNoneFactory;
import com.rtg.variant.bayes.snp.ModelSnpFactory;

import junit.framework.TestCase;

/**
 */
public class IndividualSampleFactoryTest extends TestCase {

  public void testFactory() throws Exception {
    Diagnostic.setLogStream();
    final VariantParams params = new VariantParamsBuilder()
      .genomePriors(new GenomePriorParamsBuilder().create())
      .machineErrorName("default")
      .genome(new MockReaderParams(ReaderTestUtils.getReaderDnaMemory(ReaderTestUtils.SEQ_DNA_SIMPLE)))
      .create();
    final ModelSnpFactory diploid = new ModelSnpFactory(params.genomePriors(), false, new NoAlleleBalance());
    final ModelSnpFactory haploid = new ModelSnpFactory(params.genomePriors(), true, new NoAlleleBalance());
    final ModelNoneFactory none = new ModelNoneFactory();
    final MachineErrorChooserInterface chooser = MultisampleUtils.chooser(params);
    final SexMemo sexMemo = Utils.createSexMemo(params);
    final IndividualSampleFactory<Description> individualFactories = new IndividualSampleFactory<>(params, chooser, haploid, diploid, none, params.sex(), sexMemo);
    final IndividualSampleProcessor<?> p = individualFactories.make("template", new byte[] {1, 2, 3, 4}, 2, 4);

    try {
      p.step(0);
      fail();
    } catch (IllegalStateException ex) {
      assertTrue(ex.getMessage().contains("ref=0 != base=2"));
    }

    assertNotNull(p.step(2));

  }

}
