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

package com.rtg.variant.bayes.multisample.singleton;

import java.io.File;

import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.test.BgzipFileHelper;
import com.rtg.util.test.FileHelper;
import com.rtg.variant.AlleleCountsFileConverter;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.VariantParams;
import com.rtg.variant.bayes.multisample.AbstractJointCallerConfiguration;

import htsjdk.samtools.SAMFileHeader;
import junit.framework.TestCase;

/**
 */
public class SingletonCallerConfigurationTest extends TestCase {

  public void testCreation() throws Exception {
    final File tmp = FileHelper.createTempDirectory();
    try {
      final File popFile = new File(tmp, "pop.vcf.gz");
      BgzipFileHelper.resourceToBgzipFile("com/rtg/variant/resources/pop58.vcf", popFile);

      final File alleleCountFile = new File(tmp, "allele.ac");

      final AlleleCountsFileConverter blah = new AlleleCountsFileConverter();
      blah.convert(popFile, alleleCountFile);

      Diagnostic.setLogStream();
      final VariantParams p = VariantParams.builder()
        .machineErrorName("illumina")
        .genomePriors(GenomePriorParams.builder().create())
        .populationPriors(alleleCountFile)
        .uberHeader(new SAMFileHeader())
        .create();
      final AbstractJointCallerConfiguration config = new SingletonCallerConfiguration.Configurator().getConfig(p);
      assertNotNull(config.getGenomeNames());
      assertNotNull(config.getJointCaller());
      assertEquals(1, config.numberOfGenomes());
      assertEquals(1, config.getGenomeNames().length);
      assertTrue(config.getJointCaller() instanceof SingletonCaller);

      assertNotNull(config.getOutputFormatter(p));
    } finally {
      FileHelper.deleteAll(tmp);
    }
  }

}
