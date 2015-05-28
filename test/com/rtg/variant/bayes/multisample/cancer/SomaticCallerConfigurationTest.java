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

package com.rtg.variant.bayes.multisample.cancer;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.Random;

import com.rtg.reference.Ploidy;
import com.rtg.relation.GenomeRelationships;
import com.rtg.relation.Relationship.RelationshipType;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.integrity.Exam;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantParamsBuilder;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.complex.ComplexTemplate;
import com.rtg.variant.bayes.multisample.ComplexCallerTest;
import com.rtg.variant.bayes.snp.DescriptionCommon;
import com.rtg.variant.format.VariantOutputVcfFormatter;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;

import junit.framework.TestCase;

/**
 */
public class SomaticCallerConfigurationTest extends TestCase {

  public void testCreation() throws IOException, InvalidParamsException {
    Diagnostic.setLogStream();

    final GenomeRelationships genomeRelationships = new GenomeRelationships();
    genomeRelationships.addGenome("cancer");
    genomeRelationships.addGenome("TEST");
    genomeRelationships.addRelationship(RelationshipType.ORIGINAL_DERIVED, "TEST", "cancer").setProperty("contamination", "0.3");

    final SAMFileHeader uber = ComplexCallerTest.makeHeaderWithSamples("TEST", "cancer");
    final VariantParamsBuilder b = new VariantParamsBuilder();
    b.genomePriors(GenomePriorParams.builder().create());
    b.genomeRelationships(genomeRelationships);
    b.machineErrorName("illumina");
    final VariantParams p = b.uberHeader(uber).create();

    final SomaticCallerConfiguration config = new SomaticCallerConfiguration.Configurator().getConfig(p);
    assertNotNull(config.getGenomeNames());
    assertNotNull(config.getJointCaller());
    assertEquals(2, config.numberOfGenomes());
    assertEquals(2, config.getGenomeNames().length);
    assertEquals("TEST", config.getGenomeNames()[0]);
    assertEquals("cancer", config.getGenomeNames()[1]);
    assertTrue(config.getJointCaller() instanceof ContaminatedSomaticCaller);
    assertEquals("TEST", config.getGenomeNames()[0]);
    assertEquals("cancer", config.getGenomeNames()[1]);
    assertEquals("TEST", config.getGenomeNames()[0]);
    assertEquals("cancer", config.getGenomeNames()[1]);
    assertTrue(config.handlesPloidy(Ploidy.POLYPLOID));

    final VariantOutputVcfFormatter output = config.getOutputFormatter(p);
    final ByteArrayOutputStream bos = new ByteArrayOutputStream();
    final SAMFileHeader sfh = new SAMFileHeader();
    final SAMReadGroupRecord samReadGroupRecord = new SAMReadGroupRecord("ab");
    samReadGroupRecord.setSample("cancer");
    final SAMReadGroupRecord samReadGroupRecord2 = new SAMReadGroupRecord("abc");
    samReadGroupRecord2.setSample("TEST");
    sfh.addReadGroup(samReadGroupRecord);
    sfh.addReadGroup(samReadGroupRecord2);
    output.writeHeader(bos, p, sfh);
    final String header = bos.toString();
    TestUtils.containsAll(header, "##SAMPLE=<ID=TEST,Genomes=TEST,Mixture=1.0,Description=\"Original genome\">",
        "##SAMPLE=<ID=cancer,Genomes=TEST;cancer,Mixture=0.30;0.70,Description=\"Original genome;Derived genome\">",
        "##PEDIGREE=<Derived=cancer,Original=TEST>",
        "##INFO=<ID=SOMATIC,Number=1,Type=String,Description=\"Indicates the variant is a somatic mutation\">",
        "##FORMAT=<ID=SSC,Number=1,Type=Float,Description=\"Somatic score\">");
  }


  private void checkDistribution(final double[][] a) {
    for (int i = 0; i < a.length; i++) {
      Exam.assertDistribution(a[i]);
      final double ref = a[i][i];
      for (int j = 0; j < a[i].length; j++) {
        assertTrue(a[i][j] <= ref);
      }
    }
  }

  private void gt(final double[] a, final int... k) {
    for (int i = 0; i < k.length; i += 2) {
      assertTrue(k[i] + ":" + k[i + 1], a[k[i]] > a[k[i + 1]]);
    }
  }

  private void gtCheck(final double[][] initialPriors) {
    final int a = 0;
    final int c = 1;
    final int aa = 2;
    final int cc = 3;
    final int ccc = 4;
    final int a7 = 5;
    gt(initialPriors[a], c, cc, aa, cc, aa, a7, cc, ccc);
    gt(initialPriors[c], a, aa, aa, a7, cc, aa, cc, ccc);
    gt(initialPriors[aa], a, c, cc, ccc);
    gt(initialPriors[cc], c, a, aa, a7);
    gt(initialPriors[ccc], cc, aa, cc, c);
    gt(initialPriors[a7], aa, cc, aa, a, a, c);
  }

  public void testMakeInitialPriors() {
    final byte[] template = new byte[1000];
    Arrays.fill(template, (byte) 4);
    final ComplexTemplate cot = new ComplexTemplate(template, "", 500, 502);
    final Description descr = new DescriptionCommon("A", "C", "AA", "CC", "CCC", "AAAAAAA");
    final double[][] initialPriors = SomaticCallerConfiguration.makeSomaticInitialPriors(descr, cot);
    //System.err.println(toStringLog(initialPriors));
    checkDistribution(initialPriors);
    gtCheck(initialPriors);
  }

  public void testMakeInitialPriors2() {
    final byte[] template = new byte[1000];
    Arrays.fill(template, (byte) 4);
    final ComplexTemplate cot = new ComplexTemplate(template, "", 500, 502);
    final Description descr = new DescriptionCommon("A", "AAAAAAA");
    final double[][] initialPriors = SomaticCallerConfiguration.makeSomaticInitialPriors(descr, cot);
    //System.err.println(toStringLog(initialPriors));
    checkDistribution(initialPriors);
  }

  public void testMakeInitialPriorsRandomTemplate() {
    final byte[] template = new byte[1000];
    final Random rand = new Random(42);
    for (int i = 0; i < template.length; i++) {
      template[i] = (byte) (rand.nextInt(4) + 1);
    }
    final ComplexTemplate cot = new ComplexTemplate(template, "", 500, 500);
    final Description descr = new DescriptionCommon("A", "C", "AA", "CC", "CCC", "AAAAAAA");
    final double[][] initialPriors = SomaticCallerConfiguration.makeSomaticInitialPriors(descr, cot);
    //System.err.println(toStringLog(initialPriors));
    checkDistribution(initialPriors);
    gtCheck(initialPriors);
  }
}
