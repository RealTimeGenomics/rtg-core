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

package com.rtg.variant.avr;

import java.util.Properties;

import com.rtg.vcf.VcfReader;
import com.rtg.vcf.VcfRecord;

/**
 */
public class GtQualComplexMultiplierModelBuilderTest extends AbstractModelBuilderTest<GtQualComplexMultiplierModelBuilder> {

  @Override
  GtQualComplexMultiplierModelBuilder getModelBuilder(String[] formatAttributes, String[] infoAttributes, String[] derivedAttributes) {
    return new GtQualComplexMultiplierModelBuilder(formatAttributes, infoAttributes, derivedAttributes);
  }

  public void testAnnotate() {
    final GtQualComplexMultiplierModelBuilder amb = getModelBuilder(new String[]{"GQ"}, new String[]{}, new String[]{});
    assertNotNull(amb);

    final Properties params = new Properties();
    params.setProperty(GtQualComplexMultiplierModelBuilder.PARAMETER_MULTIPLIER_COMPLEX_HETEROZYGOUS, "0.25");
    params.setProperty(GtQualComplexMultiplierModelBuilder.PARAMETER_MULTIPLIER_COMPLEX_HOMOZYGOUS, "0.5");
    params.setProperty(GtQualComplexMultiplierModelBuilder.PARAMETER_MULTIPLIER_SIMPLE_HETEROZYGOUS, "1.25");
    params.setProperty(GtQualComplexMultiplierModelBuilder.PARAMETER_MULTIPLIER_SIMPLE_HOMOZYGOUS, "1.5");
    params.setProperty(GtQualComplexMultiplierModelBuilder.PARAMETER_MULTIPLIER_QUAL_COMPLEX, "0.75");
    params.setProperty(GtQualComplexMultiplierModelBuilder.PARAMETER_MULTIPLIER_QUAL_SIMPLE, "2.75");
    amb.setModelParameters(params);

    amb.build();

    final AbstractPredictModel apm = amb.getModel();
    assertNotNull(apm);

    //System.err.println(apm.toString());

    VcfRecord record = VcfReader.vcfLineToRecord("chr5\t12041\trs55926606\tA\tT\t100\tPASS\tXRX\tGT:GQ\t0|1:12.3\t.\t1|1:99.9");
    apm.annotate(record);

    //System.err.println(record.toString());

    assertEquals(3.0, record.getSampleDouble(0, "GQ"), 0.01);
    assertEquals(50.0, record.getSampleDouble(2, "GQ"), 0.01);
    assertEquals(75.0, Double.valueOf(record.getQuality()), 0.01);


    record = VcfReader.vcfLineToRecord("chr5\t12041\trs55926606\tA\tT\t100\tPASS\t.\tGT:GQ\t0|1:12.3\t.\t1|1:99.9");
    apm.annotate(record);

    //System.err.println(record.toString());

    assertEquals(15.0, record.getSampleDouble(0, "GQ"), 0.01);
    assertEquals(150.0, record.getSampleDouble(2, "GQ"), 0.01);
    assertEquals(275.0, Double.valueOf(record.getQuality()), 0.01);
  }

  public void testDefaults() {
    final GtQualComplexMultiplierModelBuilder amb = getModelBuilder(new String[]{"GQ"}, new String[]{}, new String[] {});
    assertNotNull(amb);

    amb.build();

    final AbstractPredictModel apm = amb.getModel();
    assertNotNull(apm);

    VcfRecord record = VcfReader.vcfLineToRecord("chr5\t12041\trs55926606\tA\tT\t100\tPASS\tXRX\tGT:GQ\t0|1:12.3\t.\t1|1:99.9");
    apm.annotate(record);

    //System.err.println(record.toString());

    assertEquals(3, record.getSampleDouble(0, "GQ"), 0.01);
    assertEquals(116, record.getSampleDouble(2, "GQ"), 0.01);
    assertEquals(4, Double.valueOf(record.getQuality()), 0.01);


    record = VcfReader.vcfLineToRecord("chr5\t12041\trs55926606\tA\tT\t100\tPASS\t.\tGT:GQ\t0|1:12.3\t.\t1|1:99.9");
    apm.annotate(record);

    //System.err.println(record.toString());

    assertEquals(6, record.getSampleDouble(0, "GQ"), 0.01);
    assertEquals(200, record.getSampleDouble(2, "GQ"), 0.01);
    assertEquals(20, Double.valueOf(record.getQuality()), 0.01);

  }
}
