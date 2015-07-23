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

package com.rtg.variant.bayes.multisample.population;

import java.io.ByteArrayOutputStream;
import java.io.File;

import com.rtg.reference.Ploidy;
import com.rtg.relation.GenomeRelationships;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.test.BgzipFileHelper;
import com.rtg.util.test.FileHelper;
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
    final File tmp = FileHelper.createTempDirectory();
    try {
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
        final VariantParams params = VariantParams.builder().genomeRelationships(file).genomePriors("testhumanprior").machineErrorName("illumina").populationPriors(alleleCountFile).uberHeader(uber).create();
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
    } finally {
      FileHelper.deleteAll(tmp);
    }
  }
}
