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

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;

import com.rtg.util.TestUtils;
import com.rtg.vcf.VcfReaderTest;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.header.VcfHeader;

/**
 */
public class GtQualComplexMultiplierModelTest extends AbstractPredictModelTest<GtQualComplexMultiplierModel> {


  @Override
  GtQualComplexMultiplierModel getPredictModel() {
    return new GtQualComplexMultiplierModel();
  }

  @Override
  GtQualComplexMultiplierModel getPredictModel(InputStream is) throws IOException {
    return new GtQualComplexMultiplierModel(is);
  }

  @Override
  public void testLoadSave() throws Exception {
    final AbstractPredictModel apm = getPredictModel();
    assertNotNull(apm);
    assertEquals("AVR", apm.getField());
    assertNotNull(apm.toString());
    TestUtils.containsAll(apm.toString(),
        "multiplier.gq.simple.homozygous\t1.0",
        "multiplier.gq.simple.heterozygous\t1.0",
        "multiplier.gq.complex.homozygous\t1.0",
        "multiplier.gq.complex.heterozygous\t1.0",
        "multiplier.qual.simple\t1.0",
        "multiplier.qual.complex\t1.0");

    ((GtQualComplexMultiplierModel) apm).setGqComplexHeterozygousMultiplier(0.25);
    ((GtQualComplexMultiplierModel) apm).setGqComplexHomozygousMultiplier(0.5);
    ((GtQualComplexMultiplierModel) apm).setQualComplexMultiplier(0.75);
    ((GtQualComplexMultiplierModel) apm).setGqSimpleHeterozygousMultiplier(2.25);
    ((GtQualComplexMultiplierModel) apm).setGqSimpleHomozygousMultiplier(2.5);
    ((GtQualComplexMultiplierModel) apm).setQualSimpleMultiplier(2.75);

    TestUtils.containsAll(apm.toString(),
        "multiplier.gq.simple.homozygous\t2.5",
        "multiplier.gq.simple.heterozygous\t2.25",
        "multiplier.gq.complex.homozygous\t0.5",
        "multiplier.gq.complex.heterozygous\t0.25",
        "multiplier.qual.simple\t2.75",
        "multiplier.qual.complex\t0.75");

    try (final ByteArrayOutputStream baos = new ByteArrayOutputStream()) {
      apm.save(baos);

      //System.err.println(baos.toString());
      //System.err.println(apm.toString());

      try (final ByteArrayInputStream bais = new ByteArrayInputStream(baos.toByteArray())) {
        final AbstractPredictModel apm2 = getPredictModel(bais);
        assertNotNull(apm2);
        assertEquals("AVR", apm2.getField());
        assertEquals(apm.toString(), apm2.toString());
      }
    }
  }


  @Override
  public void testAnnotate() {
    final AbstractPredictModel apm = getPredictModel();
    assertNotNull(apm);
    assertEquals("AVR", apm.getField());

    ((GtQualComplexMultiplierModel) apm).setGqComplexHeterozygousMultiplier(0.25);
    ((GtQualComplexMultiplierModel) apm).setGqComplexHomozygousMultiplier(0.5);
    ((GtQualComplexMultiplierModel) apm).setQualComplexMultiplier(0.752);

    VcfRecord record = VcfReaderTest.vcfLineToRecord("chr5\t12041\trs55926606\tA\tT\t100\tPASS\tXRX\tGT:GQ\t0|1:12.3\t.\t1|1:99.9");
    apm.annotate(record);

    //System.err.println(record.toString());

    assertEquals(3.0, record.getSampleDouble(0, "GQ"), 0.01);
    assertEquals(50.0, record.getSampleDouble(2, "GQ"), 0.01);
    assertEquals(75.2, Double.valueOf(record.getQuality()), 0.01);

    record = VcfReaderTest.vcfLineToRecord("chr5\t12041\trs55926606\tA\tT\t100\tPASS\t.\tGT:GQ\t0|1:12.3\t.\t1|1:99.9");
    apm.annotate(record);

    //System.err.println(record.toString());

    assertEquals(12.0, record.getSampleDouble(0, "GQ"), 0.01);
    assertEquals(100.0, record.getSampleDouble(2, "GQ"), 0.01);
    assertEquals(100.0, Double.valueOf(record.getQuality()), 0.01);

    apm.annotate(record);
    apm.annotate(record);
    apm.annotate(record);


    TestUtils.containsAll(apm.getSummary(),
        "Simple homozygous\t8",
        "Simple heterozygous\t4",
        "Complex homozygous\t2",
        "Complex heterozygous\t1",
        "Total\t15",
        "Simple\t4",
        "Complex\t1",
        "Total\t5"
        );
  }

  public void testDefaultSettings() {
    final AbstractPredictModel apm = getPredictModel();
    //System.err.println(apm.toString());

    TestUtils.containsAll(apm.toString(),
        "multiplier.gq.simple.homozygous\t1.0",
        "multiplier.gq.simple.heterozygous\t1.0",
        "multiplier.gq.complex.homozygous\t1.0",
        "multiplier.gq.complex.heterozygous\t1.0",
        "multiplier.qual.simple\t1.0",
        "multiplier.qual.complex\t1.0"
        );
  }

  @Override
  public void testUpdateHeader() {
    final VcfHeader header = new VcfHeader();
    final AbstractPredictModel apm = getPredictModel();
    assertNotNull(apm);
    apm.updateHeader(header);
    //System.err.println(header);
    assertEquals(0, header.getInfoLines().size());
    assertEquals(0, header.getFormatLines().size());
  }

}
