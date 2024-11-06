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

package com.rtg.variant.bayes.multisample.statistics;

import java.io.File;

import com.rtg.launcher.MockReaderParams;
import com.rtg.launcher.OutputParams;
import com.rtg.reader.ReaderTestUtils;
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
      final OutputParams outputParams = new OutputParams(tempDir, false);
      final VariantParams p = VariantParams.builder()
        .machineErrorName("illumina")
        .genomePriors(GenomePriorParams.builder().create())
        .noComplexCalls(true)
        .outputParams(outputParams)
        .uberHeader(new SAMFileHeader())
        .genome(new MockReaderParams(ReaderTestUtils.getReaderDnaMemory(ReaderTestUtils.SEQ_DNA_SIMPLE)))
        .create();
      final AbstractJointCallerConfiguration config = new VstatsCallerConfiguration.Configurator().getConfig(p, null);
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
