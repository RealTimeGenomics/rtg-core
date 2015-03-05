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

package com.rtg.variant.bayes.multisample.family;

import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.StringReader;

import com.rtg.relation.GenomeRelationships;
import com.rtg.relation.RelationshipsFileParser;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantParamsBuilder;
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
    final VariantParamsBuilder b = new VariantParamsBuilder();
    b.genomePriors(GenomePriorParams.builder().create());
    b.genomeRelationships(rel);
    b.machineErrorName("illumina");
    final VariantParams p = b.create();

    final DiseasedFamilyCallerConfiguration config = new DiseasedFamilyCallerConfiguration.Configurator().getConfig(p, new String[] {"father", "mother", "child"});
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
