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

package com.rtg.variant.bayes.multisample.statistics;

import java.io.File;

import com.rtg.launcher.OutputParams;
import com.rtg.tabix.TabixIndexer;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.BgzipFileHelper;
import com.rtg.util.test.FileHelper;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.VariantParams;
import com.rtg.variant.bayes.multisample.AbstractJointCallerConfiguration;

import htsjdk.samtools.SAMFileHeader;
import junit.framework.TestCase;

/**
 */
public class VstatsCallerConfigurationTest extends TestCase {

  public void testCreation() throws Exception {
    final File tmp = FileHelper.createTempDirectory();
    final File tempDir = FileUtils.createTempDir("vstats", "test");
    try {
      final File popFile = new File(tmp, "pop.vcf.gz");
      BgzipFileHelper.resourceToBgzipFile("com/rtg/variant/resources/pop58.vcf", popFile);
      final TabixIndexer ti = new TabixIndexer(popFile);
      ti.saveVcfIndex();

      Diagnostic.setLogStream();
      final OutputParams outputParams = new OutputParams(tempDir, false, false);
      final VariantParams p = VariantParams.builder()
        .machineErrorName("illumina")
        .genomePriors(GenomePriorParams.builder().create())
        .noComplexCalls(true)
        .outputParams(outputParams)
        .uberHeader(new SAMFileHeader())
        .create();
      final AbstractJointCallerConfiguration config = new VstatsCallerConfiguration.Configurator().getConfig(p);
      assertNotNull(config.getGenomeNames());
      assertNotNull(config.getJointCaller());
      assertEquals(1, config.numberOfGenomes());
      assertEquals(1, config.getGenomeNames().length);
      assertTrue(config.getJointCaller() instanceof VstatsCaller);

      assertNotNull(config.getOutputFormatter(p));

      config.getJointCaller().close();

    } finally {
      assertTrue(FileHelper.deleteAll(tmp));
      assertTrue(FileHelper.deleteAll(tempDir));
    }
  }

}
