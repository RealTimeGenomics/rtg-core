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

package com.rtg.variant.bayes.multisample.population;

import java.io.ByteArrayOutputStream;
import java.io.File;

import com.rtg.launcher.MockReaderParams;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reference.Ploidy;
import com.rtg.relation.GenomeRelationships;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.BgzipFileHelper;
import com.rtg.variant.AlleleCountsFileConverter;
import com.rtg.variant.VariantParams;
import com.rtg.variant.bayes.multisample.ComplexCallerTest;
import com.rtg.variant.bayes.multisample.IndividualSampleProcessor;
import com.rtg.variant.format.VariantOutputVcfFormatter;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import junit.framework.TestCase;

/**
 */
public class PopulationCallerConfigurationTest extends TestCase {

  public void testConfigurator() throws Exception {
    try (final TestDirectory tmp = new TestDirectory()) {
      Diagnostic.setLogStream();
      final File popFile = new File(tmp, "pop.vcf.gz");
      BgzipFileHelper.resourceToBgzipFile("com/rtg/variant/resources/pop58.vcf", popFile);

      final File alleleCountFile = new File(tmp, "allele.ac");

      final AlleleCountsFileConverter blah = new AlleleCountsFileConverter();
      blah.convert(popFile, alleleCountFile);

      final MemoryPrintStream mps = new MemoryPrintStream();
      Diagnostic.setLogStream(mps.printStream());
      try {
        // Genomes are unrelated
        final GenomeRelationships file = new GenomeRelationships();
        file.addGenome("one");
        file.addGenome("two");
        file.addGenome("three");
        final SAMFileHeader uber = ComplexCallerTest.makeHeaderWithSamples("two", "one", "three");
        final VariantParams params = VariantParams.builder().genomeRelationships(file).genomePriors("testhumanprior").machineErrorName("illumina").populationPriors(alleleCountFile).uberHeader(uber).genome(new MockReaderParams(ReaderTestUtils.getReaderDnaMemory(ReaderTestUtils.SEQ_DNA_SIMPLE))).create();
        final PopulationCallerConfiguration config = new PopulationCallerConfiguration.Configurator().getConfig(params, null);
        assertNotNull(config.getOutputFormatter(params));

        assertFalse(config.getSnpHypotheses(1, null, 0).diploid().haploid());
        assertEquals(10, config.getSnpHypotheses(1, null, 0).diploid().size());

        assertTrue(config.getSnpHypotheses(1, null, 0).haploid().haploid());
        assertEquals(4, config.getSnpHypotheses(1, null, 0).haploid().size());


        assertTrue(config.getJointCaller().getClass().equals(PopulationCaller.class));

        assertTrue(config.handlesPloidy(Ploidy.POLYPLOID));
        assertTrue(config.handlesPloidy(Ploidy.HAPLOID));
        assertTrue(config.handlesPloidy(Ploidy.DIPLOID));

        final IndividualSampleProcessor<?>[] processors = config.getIndividualSampleProcessors("", new byte[] {1, 2, 3}, 0, 1);
        assertEquals(3, processors.length);

        final VariantOutputVcfFormatter output = config.getOutputFormatter(params);
        final ByteArrayOutputStream bos = new ByteArrayOutputStream();
        final SAMFileHeader sfh = new SAMFileHeader();
        final SAMReadGroupRecord samReadGroupRecord = new SAMReadGroupRecord("ab");
        samReadGroupRecord.setSample("one");
        final SAMReadGroupRecord samReadGroupRecord2 = new SAMReadGroupRecord("abc");
        samReadGroupRecord2.setSample("two");
        final SAMReadGroupRecord samReadGroupRecord3 = new SAMReadGroupRecord("abcd");
        samReadGroupRecord3.setSample("three");
        sfh.addReadGroup(samReadGroupRecord);
        sfh.addReadGroup(samReadGroupRecord2);
        sfh.addReadGroup(samReadGroupRecord3);
        output.writeHeader(bos, params, sfh);
        final String header = bos.toString();
        //System.err.println(header);
        TestUtils.containsAll(header, "one\tthree\ttwo");
      } finally {
        Diagnostic.setLogStream();
      }
      TestUtils.containsAll(mps.toString(), "Using Population caller");
    }
  }
}
