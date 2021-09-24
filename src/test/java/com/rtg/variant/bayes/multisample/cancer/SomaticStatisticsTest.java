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
