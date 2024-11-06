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

package com.rtg.variant.format;

import static com.rtg.variant.format.VcfInfoField.LOH;

import java.io.IOException;
import java.util.EnumSet;
import java.util.HashMap;

import com.rtg.calibrate.CalibratedPerSequenceExpectedCoverage;
import com.rtg.calibrate.Calibrator;
import com.rtg.calibrate.CovariateEnum;
import com.rtg.reference.Ploidy;
import com.rtg.util.TestUtils;
import com.rtg.util.intervals.ReferenceRegions;
import com.rtg.util.test.NanoRegression;
import com.rtg.variant.SomaticParamsBuilder;
import com.rtg.variant.Variant;
import com.rtg.variant.Variant.VariantFilter;
import com.rtg.variant.VariantLocus;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantSample;
import com.rtg.variant.bayes.multisample.cancer.SomaticFilterTest;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.header.VcfHeader;

import junit.framework.TestCase;

/**
 */
public class VcfInfoFieldTest extends TestCase {

  private NanoRegression mNano = null;

  @Override
  public void setUp() {
    mNano = new NanoRegression(this.getClass(), false);
  }

  @Override
  public void tearDown() throws Exception {
    try {
      mNano.finish();
    } finally {
      mNano = null;
    }
  }

  public void testEnum() {
    TestUtils.testEnum(VcfInfoField.class, "[LOH, NCS, DISEASE, RDS, DPS, DP_DPR, XRX, RCE, CT, RTRM, RSPLT, NREF, IC, EP, LAL, QD, NAA, AC, AN, SGP]");
    for (VcfInfoField field : EnumSet.range(LOH, VcfInfoField.RSPLT)) {
      assertFalse(field.isVcfAnnotator());
    }
    for (VcfInfoField field : EnumSet.range(VcfInfoField.IC, VcfInfoField.SGP)) {
      assertTrue(field.isVcfAnnotator());
    }
  }

  public void testHeaders() throws IOException {
    final VcfHeader header = new VcfHeader();
    for (VcfInfoField field : VcfInfoField.values()) {
      field.updateHeader(header);
    }
    mNano.check("vcfinfofield-header.txt", header.toString());
  }

  private final class DummyCoverageThreshold extends CalibratedPerSequenceExpectedCoverage {
    DummyCoverageThreshold() {
      super(new Calibrator(CovariateEnum.getCovariates(CovariateEnum.DEFAULT_COVARIATES, null), new ReferenceRegions()), new HashMap<>(), new HashMap<>(), null);
    }
    @Override
    public double expectedCoverage(String sequenceName, String sampleName) {
      throw new UnsupportedOperationException();
    }
    @Override
    public double expectedTotalCoverage(String sequenceName) {
      assertEquals("ref", sequenceName);
      return 0.7;
    }
  }

  public void testRecordUpdate() {
    final VariantParams params = VariantParams.builder().expectedCoverage(new DummyCoverageThreshold()).create();
    final VcfRecord record = new VcfRecord("ref", 2, "G");
    record.addAltCall("A");
    record.setNumberOfSamples(3);
    record.addFormatAndSample("GT", "1/1");
    record.addFormatAndSample("GT", "0/1");
    record.addFormatAndSample("GT", "0/0");
    final VariantSample sampleA = VariantOutputVcfFormatterTest.createSample(Ploidy.DIPLOID, "A", false, 10.0, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0);
    sampleA.setCoverage(10);
    final VariantSample sampleB = VariantOutputVcfFormatterTest.createSample(Ploidy.DIPLOID, "A:G", false, 10.0, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0);
    sampleB.setCoverage(11);
    final VariantSample sampleC = VariantOutputVcfFormatterTest.createSample(Ploidy.DIPLOID, "G", true, 10.0, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0);
    sampleC.setCoverage(12);
    final Variant call = new Variant(new VariantLocus("ref", 2, 3, "G", 'C'), sampleA, sampleB, sampleC);
    call.setPossibleCause("A");
    call.setPossibleCauseScore(15.0);
    call.setNormalCancerScore(7.5);
    call.setDiseasePresenceScore(8.5);
    call.setComplexEquivalent();
    call.setComplexScored();
    call.addFilter(VariantFilter.COVERAGE);
    call.setTrimmed();
    call.setSplitId(1);
    VcfFormatField.DP.updateRecordSample(record, call, sampleA, null, params, false);
    VcfFormatField.DP.updateRecordSample(record, call, sampleB, null, params, false);
    VcfFormatField.DP.updateRecordSample(record, call, sampleC, null, params, false);
    record.setQuality("123.4");
    for (VcfInfoField field : VcfInfoField.values()) {
      field.updateRecord(record, call, params, false);
    }
    assertEquals("ref\t3\t.\tG\tA\t123.4\t.\tNCS=32.574;DISEASE=A;RDS=6.5;DPS=36.9;DP=33;DPR=47.143;XRX;RCE;CT=2147483647;RTRM;RSPLT=1;IC=0.333;EP=0.724;LAL=1;QD=3.739;NAA=1;AC=3;AN=6\tGT:DP\t1/1:10\t0/1:11\t0/0:12", record.toString());
  }

  public void testMultiAlleleAC() {
    final VcfRecord record = new VcfRecord("ref", 2, "G");
    record.setNumberOfSamples(3);
    record.addAltCall("A");
    record.addAltCall("C");
    record.addFormatAndSample("GT", "1/1");
    record.addFormatAndSample("GT", "2/1");
    record.addFormatAndSample("GT", "0/2");
    VcfInfoField.AC.updateRecord(record, null, null, false);
    assertEquals("ref\t3\t.\tG\tA,C\t.\t.\tAC=3,2\tGT\t1/1\t2/1\t0/2", record.toString());
  }

  public void testNoCoverage() {
    final VcfRecord record = new VcfRecord("ref", 2, "A");
    VcfInfoField.DP_DPR.updateRecord(record, new Variant(new VariantLocus("ref", 2, 3, "A", 'G')), VariantParams.builder().create(), false);
    assertEquals("ref\t3\t.\tA\t.\t.\t.\t.", record.toString());
  }

  public void testDummyCoverageThreshold() {
    final VcfRecord record = new VcfRecord("ref", 2, "A");
    final Variant variant = new Variant(new VariantLocus("ref", 2, 3, "A", 'G'));
    variant.addFilter(VariantFilter.COVERAGE);
    VcfInfoField.CT.updateRecord(record, variant, null, false);
    assertEquals("ref\t3\t.\tA\t.\t.\t.\tCT=0", record.toString());
  }

  public void testFormatPossibleCause() {
    final Variant call = new Variant(new VariantLocus("ref", 2, 3, "A", 'C'));
    assertEquals("*", VcfInfoField.formatPossibleCause(call, false));
    call.setPossibleCause("G");
    assertEquals("G", VcfInfoField.formatPossibleCause(call, false));
    assertEquals("CG", VcfInfoField.formatPossibleCause(call, true));
    call.setPossibleCause("G:T");
    assertEquals("CG:CT", VcfInfoField.formatPossibleCause(call, true));
  }
  public void testLohInfo() {
    final VcfRecord record = SomaticFilterTest.getSomaticVcfRecord("0/1", "0/0");
    final VariantParams variantParams = VariantParams.builder().somaticParams(new SomaticParamsBuilder().lohPrior(0.01).create()).create();
    LOH.updateRecord(record, null, variantParams,  false);
    assertEquals("1.0", record.getInfo("LOH"));
  }
  public void testLohNotEnabled() {
    final VcfRecord record = SomaticFilterTest.getSomaticVcfRecord("0/1", "0/0");
    final VariantParams variantParams = VariantParams.builder().somaticParams(new SomaticParamsBuilder().lohPrior(0.0).create()).create();
    LOH.updateRecord(record, null, variantParams,  false);
    assertNull(record.getInfo("LOH"));
  }
}
