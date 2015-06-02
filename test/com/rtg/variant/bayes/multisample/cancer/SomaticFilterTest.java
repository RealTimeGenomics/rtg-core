/*
 * Copyright (c) 2015. Real Time Genomics Limited.
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

import com.rtg.relation.GenomeRelationships;
import com.rtg.relation.Relationship;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.VariantParams;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.header.VcfHeader;

import junit.framework.TestCase;

/**
 * Tests the corresponding class.
 */
public class SomaticFilterTest extends TestCase {

  public void test() {
    final GenomeRelationships genomeRelationships = new GenomeRelationships();
    genomeRelationships.addGenome("cancer");
    genomeRelationships.addGenome("normal");
    genomeRelationships.addRelationship(Relationship.RelationshipType.ORIGINAL_DERIVED, "normal", "cancer").setProperty("contamination", "0.3");
    final VariantParams params = VariantParams.builder()
      .genomePriors(GenomePriorParams.builder().create())
      .genomeRelationships(genomeRelationships)
      .create();
    final SomaticStatistics ss = new SomaticStatistics(params, "GQ");
    final SomaticFilter sf = new SomaticFilter(ss, false);
    final VcfHeader vcfHeader = new VcfHeader();
    vcfHeader.addSampleName("cancer");
    vcfHeader.addSampleName("normal");
    sf.setHeader(vcfHeader);
    final VcfRecord rec1 = new VcfRecord();
    rec1.setRefCall("A")
      .setNumberOfSamples(2)
      .setSequence("pretend")
      .setStart(42)
      .addFormatAndSample("GT", "0/0")
      .addFormatAndSample("GT", "0/1")
      .addFormatAndSample("GQ", "20")
      .addFormatAndSample("GQ", "20")
      .addFormatAndSample("AD", "50,0")
      .addFormatAndSample("AD", "25,25");
    assertFalse(sf.accept(rec1));
    rec1.setInfo("SOMATIC", "A:C");
    assertTrue(sf.accept(rec1));
  }
}
