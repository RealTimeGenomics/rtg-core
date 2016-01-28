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

package com.rtg.variant.bayes.multisample;


import com.rtg.launcher.MockReaderParams;
import com.rtg.mode.SequenceMode;
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
      .genome(new MockReaderParams(ReaderTestUtils.getReaderDnaMemory(ReaderTestUtils.SEQ_DNA_SIMPLE), SequenceMode.UNIDIRECTIONAL))
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
