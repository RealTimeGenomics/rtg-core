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

  public static VcfRecord getSomaticVcfRecord(String normalGt, String cancerGt) {
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

  static VcfRecord getVcfRecord(String normalGt, String cancerGt) {
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
