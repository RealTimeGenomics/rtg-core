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
import java.io.StringReader;

import com.rtg.launcher.MockReaderParams;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.relation.GenomeRelationships;
import com.rtg.relation.RelationshipsFileParser;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
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
public class DiseasedFamilyCallerConfigurationTest extends TestCase {

  public void testCreation() throws Exception {
    Diagnostic.setLogStream();
    final GenomeRelationships rel = RelationshipsFileParser.load(new BufferedReader(new StringReader(DiseasedFamilyPosteriorTest.RELATIONS1)));
    final SAMFileHeader uber = ComplexCallerTest.makeHeaderWithSamples("father", "mother", "child");
    final VariantParamsBuilder b = new VariantParamsBuilder();
    b.genomePriors(GenomePriorParams.builder().create());
    b.genomeRelationships(rel);
    b.machineErrorName("illumina");
    b.uberHeader(uber);
    b.genome(new MockReaderParams(ReaderTestUtils.getReaderDnaMemory(ReaderTestUtils.SEQ_DNA_SIMPLE)));
    final VariantParams p = b.create();

    final DiseasedFamilyCallerConfiguration config = new DiseasedFamilyCallerConfiguration.Configurator().getConfig(p, null);
    assertNotNull(config.getGenomeNames());
    assertNotNull(config.getJointCaller());
    assertEquals(3, config.numberOfGenomes());
    assertEquals(3, config.getGenomeNames().length);
    assertEquals("father", config.getGenomeNames()[0]);
    assertEquals("mother", config.getGenomeNames()[1]);
    assertEquals("child", config.getGenomeNames()[2]);
    assertTrue(config.getJointCaller() instanceof DiseasedFamilyCaller);

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
    TestUtils.containsAll(header, "##diseased=mother",
        "##PEDIGREE=<Child=child,Mother=mother,Father=father>",
        "##INFO=<ID=DISEASE,Number=1,Type=String,Description=\"Indicates the variant is linked to the disease\">",
        "##INFO=<ID=RDS,Number=1,Type=Float,Description=\"RTG disease call score\">");
  }
}
