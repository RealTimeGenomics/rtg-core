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

package com.rtg.variant.bayes.multisample.family;

import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.StringReader;

import com.rtg.launcher.MockReaderParams;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reference.Ploidy;
import com.rtg.relation.GenomeRelationships;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.test.BgzipFileHelper;
import com.rtg.util.test.FileHelper;
import com.rtg.variant.AlleleCountsFileConverter;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantParamsBuilder;
import com.rtg.variant.bayes.multisample.ComplexCallerTest;
import com.rtg.variant.format.VariantOutputVcfFormatter;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import junit.framework.TestCase;

/**
 */
public class FamilyCallerConfigurationTest extends TestCase {

  static final String TRIO_PED = ""
      + "0\tfather\t0\t0\t1\t0" + StringUtils.LS
      + "0\tmother\t0\t0\t2\t0" + StringUtils.LS
      + "0\tchild\tfather\tmother\t1\t0" + StringUtils.LS;

  static final String TOOMANY_PED = TRIO_PED
      + "0\tchild2\tfather\tanothermother\t1\t0" + StringUtils.LS;


  public void testCreation() throws Exception {
    final File tmp = FileHelper.createTempDirectory();
    try {
      final File popFile = new File(tmp, "pop.vcf.gz");
      BgzipFileHelper.resourceToBgzipFile("com/rtg/variant/resources/pop58.vcf", popFile);

      final File alleleCountFile = new File(tmp, "allele.ac");

      final AlleleCountsFileConverter blah = new AlleleCountsFileConverter();
      blah.convert(popFile, alleleCountFile);

      Diagnostic.setLogStream();

      try {
        final GenomeRelationships rel = GenomeRelationships.loadGenomeRelationships(new BufferedReader(new StringReader(TOOMANY_PED)));
        final VariantParams p = new VariantParamsBuilder().genomeRelationships(rel).uberHeader(ComplexCallerTest.makeHeaderWithSamples("child", "mother", "father")).create();
        new FamilyCallerConfiguration.Configurator().getConfig(p, null);
        fail("Should not have liked the family");
      } catch (NoTalkbackSlimException e) {
        // Expected
      }

      final GenomeRelationships rel = GenomeRelationships.loadGenomeRelationships(new BufferedReader(new StringReader(TRIO_PED)));
      final VariantParamsBuilder b = new VariantParamsBuilder();
      b.genomePriors(GenomePriorParams.builder().create());
      b.genomeRelationships(rel);
      b.machineErrorName("illumina");
      b.populationPriors(alleleCountFile);
      b.uberHeader(ComplexCallerTest.makeHeaderWithSamples("child", "mother", "father"));
      b.genome(new MockReaderParams(ReaderTestUtils.getReaderDnaMemory(ReaderTestUtils.SEQ_DNA_SIMPLE)));
      final VariantParams p = b.create();

      final FamilyCallerConfiguration config = new FamilyCallerConfiguration.Configurator().getConfig(p, null);
      assertNotNull(config.getGenomeNames());
      assertNotNull(config.getJointCaller());
      assertEquals(3, config.numberOfGenomes());
      assertEquals(3, config.getGenomeNames().length);
      assertEquals("father", config.getGenomeNames()[0]);
      assertEquals("mother", config.getGenomeNames()[1]);
      assertEquals("child", config.getGenomeNames()[2]);
      assertTrue(config.getJointCaller() instanceof FamilyCaller);
      assertTrue(config.handlesPloidy(Ploidy.POLYPLOID));
      assertTrue(config.handlesPloidy(Ploidy.HAPLOID));
      assertTrue(config.handlesPloidy(Ploidy.DIPLOID));

      final VariantOutputVcfFormatter output = config.getOutputFormatter(p);
      final ByteArrayOutputStream bos = new ByteArrayOutputStream();
      final SAMFileHeader sfh = new SAMFileHeader();
      final SAMReadGroupRecord samReadGroupRecord = new SAMReadGroupRecord("ab");
      samReadGroupRecord.setSample("father");
      final SAMReadGroupRecord samReadGroupRecord2 = new SAMReadGroupRecord("abc");
      samReadGroupRecord2.setSample("mother");
      final SAMReadGroupRecord samReadGroupRecord3 = new SAMReadGroupRecord("abcd");
      samReadGroupRecord3.setSample("child");
      sfh.addReadGroup(samReadGroupRecord);
      sfh.addReadGroup(samReadGroupRecord2);
      sfh.addReadGroup(samReadGroupRecord3);
      output.writeHeader(bos, p, sfh);
      final String header = bos.toString();
      //System.err.println(header);
      TestUtils.containsAll(header, "##PEDIGREE=<Child=child,Mother=mother,Father=father>");
    } finally {
      FileHelper.deleteAll(tmp);
    }
  }

  public void testBadRelationships() throws Exception {
    Diagnostic.setLogStream();

    final GenomeRelationships rel = GenomeRelationships.loadGenomeRelationships(new BufferedReader(new StringReader(TRIO_PED)));
    final VariantParamsBuilder b = new VariantParamsBuilder();
    b.genomePriors(GenomePriorParams.builder().create());
    b.genomeRelationships(rel);
    b.machineErrorName("illumina");
    b.uberHeader(ComplexCallerTest.makeHeaderWithSamples("child", "mother", "father"));
    b.genome(new MockReaderParams(ReaderTestUtils.getReaderDnaMemory(ReaderTestUtils.SEQ_DNA_SIMPLE)));
    final VariantParams p = b.create();

    FamilyCallerConfiguration config = new FamilyCallerConfiguration.Configurator().getConfig(p, null);
    assertEquals(3, config.numberOfGenomes());

    config = new FamilyCallerConfiguration.Configurator().getConfig(b.uberHeader(ComplexCallerTest.makeHeaderWithSamples("child", "mother")).create(), null);
    assertEquals(3, config.numberOfGenomes());

    config = new FamilyCallerConfiguration.Configurator().getConfig(b.uberHeader(ComplexCallerTest.makeHeaderWithSamples("child", "father")).create(), null);
    assertEquals(3, config.numberOfGenomes());

    try {
      new FamilyCallerConfiguration.Configurator().getConfig(b.uberHeader(ComplexCallerTest.makeHeaderWithSamples("father", "mother")).create(), null);
      fail("Accepted childless family");
    } catch (NoTalkbackSlimException e) {
      assertEquals("Not enough family members have mapping data provided", e.getMessage());
    }

    try {
      new FamilyCallerConfiguration.Configurator().getConfig(b.uberHeader(ComplexCallerTest.makeHeaderWithSamples("father")).create(), null);
      fail("Accepted father only family");
    } catch (NoTalkbackSlimException e) {
      assertEquals("Not enough family members have mapping data provided", e.getMessage());
    }

    try {
      new FamilyCallerConfiguration.Configurator().getConfig(b.uberHeader(ComplexCallerTest.makeHeaderWithSamples("child")).create(), null);
      fail("Accepted child only family");
    } catch (NoTalkbackSlimException e) {
      assertEquals("Not enough family members have mapping data provided", e.getMessage());
    }
  }
}
