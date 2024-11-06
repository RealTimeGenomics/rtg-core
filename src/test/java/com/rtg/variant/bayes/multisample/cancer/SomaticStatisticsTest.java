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

import java.util.Arrays;

import com.rtg.reference.Ploidy;
import com.rtg.relation.GenomeRelationships;
import com.rtg.relation.Relationship.RelationshipType;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.Variant;
import com.rtg.variant.VariantLocus;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantParamsBuilder;
import com.rtg.variant.VariantSample;
import com.rtg.variant.avr.AbstractPredictModel;
import com.rtg.variant.bayes.MockGenotypeMeasure;
import com.rtg.variant.bayes.NoNonIdentityMeasure;
import com.rtg.variant.format.VariantOutputVcfFormatter;
import com.rtg.vcf.VcfRecord;

import junit.framework.TestCase;

/**
 */
public class SomaticStatisticsTest extends TestCase {

  public void testNoddy() throws Exception {
    final GenomeRelationships genomeRelationships = new GenomeRelationships();
    genomeRelationships.addGenome("cancer");
    genomeRelationships.addGenome("TEST");
    genomeRelationships.addRelationship(RelationshipType.ORIGINAL_DERIVED, "TEST", "cancer").setProperty("contamination", "0.3");

    final VariantParamsBuilder b = new VariantParamsBuilder();
    b.genomePriors(GenomePriorParams.builder().create());
    b.genomeRelationships(genomeRelationships);
    b.machineErrorName("illumina");
    final VariantParams p = b.create();

    final MemoryPrintStream mps = new MemoryPrintStream();
    final SomaticStatistics mss = new SomaticStatistics(p, AbstractPredictModel.AVR);
    mss.printStatistics(mps.outputStream());
    assertEquals(""
        + "Passed Filters               : 0" + StringUtils.LS
        + "Failed Filters               : 0" + StringUtils.LS, mps.toString());

    mps.reset();

    final VariantSample vs1 = new VariantSample(Ploidy.DIPLOID, "G:G", false, new NoNonIdentityMeasure(new MockGenotypeMeasure(3.0)), VariantSample.DeNovoStatus.UNSPECIFIED, null);
    vs1.setStatisticsString("\tG\t4\t0.020");
    vs1.setCoverage(4, 0.020);
    final VariantSample vs2 = new VariantSample(Ploidy.DIPLOID, "T:G", false, new NoNonIdentityMeasure(new MockGenotypeMeasure(3.0)), VariantSample.DeNovoStatus.UNSPECIFIED, null);
    vs2.setStatisticsString("\tT\t4\t0.020\tG\t4\t0.020");
    vs2.setCoverage(4, 0.020);

    final VariantLocus locus = new VariantLocus("test", 4, 5, "A", 'N');
    final Variant v = new Variant(locus, vs1, vs2);

    final String[] sampleNames = {"TEST", "cancer"};
    final VcfRecord rec = new VariantOutputVcfFormatter(sampleNames).makeVcfRecord(v);

    mss.tallyVariant(rec, Arrays.asList(sampleNames));
    TestUtils.containsAll(mss.getStatistics(),
            "Passed Filters               : 1",
            "Failed Filters               : 0",
            "SNPs                         : 1",
            "MNPs                         : 0",
            "Insertions                   : 0",
            "Deletions                    : 0",
            "Indels                       : 0",
            "SNP Transitions/Transversions: 1.00 (1/1)",
            "Total Het/Hom ratio          : - (1/0)",
            "SNP Het/Hom ratio            : - (1/0)",
            "MNP Het/Hom ratio            : - (0/0)",
            "Insertion Het/Hom ratio      : - (0/0)",
            "Deletion Het/Hom ratio       : - (0/0)",
            "Indel Het/Hom ratio          : - (0/0)",
            "Insertion/Deletion ratio     : - (0/0)",
            "Indel/SNP+MNP ratio          : 0.00 (0/1)");

  }
}
