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
    final SomaticFilter sf = getSomaticFilter(true, true, true);
    final VcfRecord record = getVcfRecord("0/0", "0/1");
    assertFalse(sf.accept(record));
    record.addFormatAndSample("SS", "0")
      .addFormatAndSample("SS", "2");
    assertTrue(sf.accept(record));
  }

  public void testLossOfHeterozygosity() {
    final SomaticFilter sf = getSomaticFilter(true, false, true);
    final VcfRecord record = getSomaticVcfRecord("0/1", "0/0");
    assertFalse(sf.accept(record));
  }

  public void testLossOfRefHeterozygosity() {
    final SomaticFilter sf = getSomaticFilter(true, false, true);
    final VcfRecord record = getSomaticVcfRecord("0/1", "1/1");
    assertFalse(sf.accept(record));
  }


  public void testExcludeGainOfRefRejected() {
    final SomaticFilter sf = getSomaticFilter(true, true, false);
    final VcfRecord record = getSomaticVcfRecord("0/1", "0/0");
    assertFalse(sf.accept(record));
  }

  public void testExcludeGainOfRefAccepted() {
    final SomaticFilter sf = getSomaticFilter(true, true, false);
    final VcfRecord record = getSomaticVcfRecord("0/1", "1/1");
    assertTrue(sf.accept(record));
  }

  private VcfRecord getSomaticVcfRecord(String normalGt, String cancerGt) {
    final VcfRecord record = getVcfRecord(normalGt, cancerGt);
    record.addFormatAndSample("SS", "0")
      .addFormatAndSample("SS", "2");
    return record;
  }


  private SomaticFilter getSomaticFilter(boolean somaticOnly, boolean lossOfHeterozygosity, boolean gainOfReference) {
    final GenomeRelationships genomeRelationships = new GenomeRelationships();
    genomeRelationships.addGenome("cancer");
    genomeRelationships.addGenome("normal");
    genomeRelationships.addRelationship(Relationship.RelationshipType.ORIGINAL_DERIVED, "normal", "cancer").setProperty("contamination", "0.3");
    final VariantParams params = VariantParams.builder()
      .genomePriors(GenomePriorParams.builder().create())
      .genomeRelationships(genomeRelationships)
      .create();
    final SomaticStatistics ss = new SomaticStatistics(params, "GQ");
    final SomaticFilter sf = new SomaticFilter(ss, somaticOnly, lossOfHeterozygosity, gainOfReference);
    final VcfHeader vcfHeader = new VcfHeader();
    vcfHeader.addSampleName("cancer");
    vcfHeader.addSampleName("normal");
    sf.setHeader(vcfHeader);
    return sf;
  }

  private VcfRecord getVcfRecord(String normalGt, String cancerGt) {
    final VcfRecord record = new VcfRecord("pretend", 42, "A");
    record.setNumberOfSamples(2)
      .addFormatAndSample("GT", normalGt)
      .addFormatAndSample("GT", cancerGt)
      .addFormatAndSample("GQ", "20")
      .addFormatAndSample("GQ", "20")
      .addFormatAndSample("AD", "50,0")
      .addFormatAndSample("AD", "25,25");
    return record;
  }
}
