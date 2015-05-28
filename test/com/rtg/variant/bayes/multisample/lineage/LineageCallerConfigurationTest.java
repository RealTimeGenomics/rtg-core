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
package com.rtg.variant.bayes.multisample.lineage;

import java.io.ByteArrayOutputStream;

import com.rtg.reference.Ploidy;
import com.rtg.relation.GenomeRelationships;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.variant.VariantParams;
import com.rtg.variant.bayes.multisample.ComplexCallerTest;
import com.rtg.variant.bayes.multisample.IndividualSampleProcessor;
import com.rtg.variant.bayes.multisample.MultisampleJointCaller;
import com.rtg.variant.format.VariantOutputVcfFormatter;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;

import junit.framework.TestCase;

/**
 */
public class LineageCallerConfigurationTest extends TestCase {
  public void testConfigurator() throws Exception {
      final MemoryPrintStream mps = new MemoryPrintStream();
      Diagnostic.setLogStream(mps.printStream());
      try {
        // Genomes are unrelated
        final GenomeRelationships file = new GenomeRelationships();
        file.addGenome("zero");
        file.addGenome("one");
        file.addGenome("two");
        file.addGenome("three");
        file.addParentChild("zero", "one");
        file.addParentChild("one", "two");
        file.addParentChild("zero", "three");
        final SAMFileHeader uber = ComplexCallerTest.makeHeaderWithSamples("zero", "one", "two", "three");
        final VariantParams params = VariantParams.builder().genomeRelationships(file).genomePriors("testhumanprior").machineErrorName("illumina").uberHeader(uber).create();
        final LineageCallerConfiguration config = new LineageCallerConfiguration.Configurator().getConfig(params);
        assertNotNull(config.getOutputFormatter(params));

        assertFalse(config.getSnpHypotheses(1, null, 0).diploid().haploid());
        assertEquals(10, config.getSnpHypotheses(1, null, 0).diploid().size());

        assertTrue(config.getSnpHypotheses(1, null, 0).haploid().haploid());
        assertEquals(4, config.getSnpHypotheses(1, null, 0).haploid().size());


        final MultisampleJointCaller jointCaller = config.getJointCaller();
        assertTrue(jointCaller.getClass().equals(Lineage.class));
        final Lineage lineage = (Lineage) jointCaller;
        assertTrue(lineage.isRoot(3));
        // int values dictated by out of order used in config parameter above
        assertEquals(3, lineage.parent(0));
        assertEquals(0, lineage.parent(2));
        assertEquals(3, lineage.parent(1));

        assertFalse(config.handlesPloidy(Ploidy.POLYPLOID));
        assertTrue(config.handlesPloidy(Ploidy.HAPLOID));
        assertTrue(config.handlesPloidy(Ploidy.DIPLOID));

        final IndividualSampleProcessor<?>[] processors = config.getIndividualSampleProcessors("", new byte[] {1, 2, 3}, 0, 1);
        assertEquals(4, processors.length);

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
        TestUtils.containsAll(header, "one\tthree\ttwo\tzero");
      } finally {
        Diagnostic.setLogStream();
      }
      TestUtils.containsAll(mps.toString(), "Using Lineage caller");
  }
  public void testMultiparentFail() throws Exception {
    final MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());
    try {
      // Genomes are unrelated
      final GenomeRelationships file = new GenomeRelationships();
      file.addGenome("zero");
      file.addGenome("one");
      file.addGenome("two");
      file.addParentChild("zero", "two");
      file.addParentChild("one", "two");
      final SAMFileHeader uber = ComplexCallerTest.makeHeaderWithSamples("zero", "two", "one");
      final VariantParams params = VariantParams.builder().genomeRelationships(file).genomePriors("testhumanprior").machineErrorName("illumina").uberHeader(uber).create();
      try {
        new LineageCallerConfiguration.Configurator().getConfig(params);
        fail();
      } catch (NoTalkbackSlimException e) {
        //expected
      }

    } finally {
      Diagnostic.setLogStream();
    }
    TestUtils.containsAll(mps.toString(), "Using Lineage caller");
  }
  public void testNoRelationshipFail() throws Exception {
    final MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());
    try {
      // Genomes are unrelated
      final GenomeRelationships file = new GenomeRelationships();
      file.addGenome("zero");
      file.addGenome("one");
      file.addGenome("two");
      final SAMFileHeader uber = ComplexCallerTest.makeHeaderWithSamples("zero", "two", "one");
      final VariantParams params = VariantParams.builder().genomeRelationships(file).genomePriors("testhumanprior").machineErrorName("illumina").uberHeader(uber).create();
      try {
        new LineageCallerConfiguration.Configurator().getConfig(params);
        fail();
      } catch (NoTalkbackSlimException e) {
        //expected
      }

    } finally {
      Diagnostic.setLogStream();
    }
    TestUtils.containsAll(mps.toString(), "Using Lineage caller");
  }
}
