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

package com.rtg.variant.avr;

import java.util.Properties;

import com.rtg.vcf.VcfReaderTest;
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

    VcfRecord record = VcfReaderTest.vcfLineToRecord("chr5\t12041\trs55926606\tA\tT\t100\tPASS\tXRX\tGT:GQ\t0|1:12.3\t.\t1|1:99.9");
    apm.annotate(record);

    //System.err.println(record.toString());

    assertEquals(3.0, record.getSampleDouble(0, "GQ"), 0.01);
    assertEquals(50.0, record.getSampleDouble(2, "GQ"), 0.01);
    assertEquals(75.0, Double.valueOf(record.getQuality()), 0.01);


    record = VcfReaderTest.vcfLineToRecord("chr5\t12041\trs55926606\tA\tT\t100\tPASS\t.\tGT:GQ\t0|1:12.3\t.\t1|1:99.9");
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

    VcfRecord record = VcfReaderTest.vcfLineToRecord("chr5\t12041\trs55926606\tA\tT\t100\tPASS\tXRX\tGT:GQ\t0|1:12.3\t.\t1|1:99.9");
    apm.annotate(record);

    //System.err.println(record.toString());

    assertEquals(3, record.getSampleDouble(0, "GQ"), 0.01);
    assertEquals(116, record.getSampleDouble(2, "GQ"), 0.01);
    assertEquals(4, Double.valueOf(record.getQuality()), 0.01);


    record = VcfReaderTest.vcfLineToRecord("chr5\t12041\trs55926606\tA\tT\t100\tPASS\t.\tGT:GQ\t0|1:12.3\t.\t1|1:99.9");
    apm.annotate(record);

    //System.err.println(record.toString());

    assertEquals(6, record.getSampleDouble(0, "GQ"), 0.01);
    assertEquals(200, record.getSampleDouble(2, "GQ"), 0.01);
    assertEquals(20, Double.valueOf(record.getQuality()), 0.01);

  }
}
