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

package com.rtg.variant.format;

import static com.rtg.variant.format.VcfFormatField.AQ;
import static com.rtg.variant.format.VcfFormatField.QA;
import static com.rtg.variant.format.VcfFormatField.VA;

import java.io.IOException;
import java.util.Collections;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import com.rtg.calibrate.CalibratedPerSequenceExpectedCoverage;
import com.rtg.calibrate.Calibrator;
import com.rtg.calibrate.Covariate;
import com.rtg.calibrate.CovariateEnum;
import com.rtg.calibrate.CovariateReadGroup;
import com.rtg.calibrate.CovariateSequence;
import com.rtg.launcher.AbstractNanoTest;
import com.rtg.reference.Ploidy;
import com.rtg.util.TestUtils;
import com.rtg.util.intervals.ReferenceRegions;
import com.rtg.util.intervals.RegionRestriction;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.Variant;
import com.rtg.variant.Variant.VariantFilter;
import com.rtg.variant.VariantLocus;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantParamsBuilder;
import com.rtg.variant.VariantSample;
import com.rtg.variant.bayes.AlleleStatisticsInt;
import com.rtg.variant.bayes.ArrayGenotypeMeasure;
import com.rtg.variant.bayes.GenotypeMeasure;
import com.rtg.variant.bayes.MockGenotypeMeasure;
import com.rtg.variant.bayes.snp.DescriptionSnp;
import com.rtg.variant.bayes.snp.HypothesesSnp;
import com.rtg.variant.bayes.snp.StatisticsSnp;
import com.rtg.variant.util.arithmetic.SimplePossibility;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.header.VcfHeader;

/**
 */
public class VcfFormatFieldTest extends AbstractNanoTest {

  public void testEnum() {
    TestUtils.testEnum(VcfFormatField.class, "[GT, VA, DP, DPR, RE, AR, RQ, GQ, RP, DN, DNP, ABP, SBP, RPB, PPB, AQ, PUR, RS, ADE, AD, SSC, SS, GL, GQD, ZY, PD, SCONT, VAF, QA, MEANQAD]");
    for (VcfFormatField field : EnumSet.range(VcfFormatField.GT, VcfFormatField.AD)) {
      assertFalse(field.isVcfAnnotator());
    }
    for (VcfFormatField field : EnumSet.range(VcfFormatField.GQD, VcfFormatField.SCONT)) {
      assertTrue(field.isVcfAnnotator());
    }
    assertEquals(0, VcfFormatField.GT.ordinal());
  }

  public void testHeaders() throws IOException {
    final VcfHeader header = new VcfHeader();
    for (VcfFormatField field : VcfFormatField.values()) {
      field.updateHeader(header);
    }
    mNano.check("vcfformatfield-headers.vcf", header.toString());
  }

  private final class DummyCoverageThreshold extends CalibratedPerSequenceExpectedCoverage {
    DummyCoverageThreshold() {
      super(new Calibrator(CovariateEnum.getCovariates(CovariateEnum.DEFAULT_COVARIATES, null), new ReferenceRegions()), new HashMap<String, Integer>(), new HashMap<String, String>(), null);
    }
    @Override
    public double expectedCoverage(String sequenceName, String sampleName) {
      assertEquals("ref", sequenceName);
      assertEquals("SAMPLE", sampleName);
      return 0.7;
    }
    @Override
    public double expectedTotalCoverage(String sequenceName) {
      throw new UnsupportedOperationException();
    }
  }

  public void testHasValue() {
    final String[] sampleNames = {"SAMPLE"};
    Variant call = new Variant(new VariantLocus("A", 2, 3, "A", 'C'));
    final VcfRecord rec = new VcfRecord("A", 2, "A");
    rec.setNumberOfSamples(sampleNames.length);
    final VariantParams params = VariantParams.builder().expectedCoverage(new DummyCoverageThreshold()).create();
    assertTrue(VcfFormatField.GT.hasValue(rec, call, null, null, params));
    for (VcfFormatField field : EnumSet.range(VcfFormatField.DP, VcfFormatField.PD)) {
      assertFalse(field.hasValue(rec, call, null, null, params));
    }
    VariantSample sample = new VariantSample(Ploidy.NONE);
    for (VcfFormatField field : EnumSet.range(VcfFormatField.DP, VcfFormatField.PD)) {
      assertFalse(field.hasValue(rec, call, sample, null, params));
    }
    final double[] measures = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.1, 0.1};
    final GenotypeMeasure measure = new ArrayGenotypeMeasure(SimplePossibility.SINGLETON, measures, new HypothesesSnp(SimplePossibility.SINGLETON, GenomePriorParams.builder().create(), false, 0));
    sample = new VariantSample(Ploidy.DIPLOID, "A:G", false, measure, VariantSample.DeNovoStatus.IS_DE_NOVO, 10.0);
    sample.setCoverage(10);
    sample.setCoverageCorrection(5.5);
    sample.setHoeffdingAlleleBalanceHet(1.0);
    sample.setHoeffdingAlleleBalanceHom(2.0);
    sample.setHoeffdingStrandBiasAllele1(3.0);
    sample.setHoeffdingStrandBiasAllele2(4.0);
    sample.setHoeffdingUnmatedBiasAllele1(5.0);
    sample.setHoeffdingUnmatedBiasAllele2(6.0);
    sample.setHoeffdingReadPositionBias(7.0);
    sample.setAmbiguityRatio(5.0);
    sample.setStatisticsString("Q");
    sample.setStats(new StatisticsSnp(DescriptionSnp.SINGLETON));
    sample.setPlacedUnmappedRatio(0.42);
    sample.setVariantAllele("C");

    call = new Variant(new VariantLocus("ref", 2, 3, "A", 'C'), sample);
    for (VcfFormatField field : EnumSet.range(VcfFormatField.GT, VcfFormatField.AD)) {
      assertTrue(field.name(), field.hasValue(rec, call, sample, sampleNames[0], params));
      field.updateRecord(rec, call, sampleNames, params, false);
    }
    for (VcfFormatField field : EnumSet.range(VcfFormatField.GQD, VcfFormatField.PD)) {
      assertTrue(field.name(), field.hasValue(rec, call, sample, sampleNames[0], params));
    }
    assertFalse(VcfFormatField.DPR.hasValue(rec, call, sample, null, params));
    assertFalse(VcfFormatField.DPR.hasValue(rec, call, sample, sampleNames[0], null));
    assertFalse(VcfFormatField.DPR.hasValue(rec, call, sample, sampleNames[0], VariantParams.builder().create()));
    call.setComplexScored();
    assertTrue(VcfFormatField.AR.hasValue(rec, call, sample, null, params));
    assertFalse(VcfFormatField.RS.hasValue(rec, call, sample, null, params));
    call.addFilter(VariantFilter.FAILED_COMPLEX);
    assertFalse(VcfFormatField.DN.hasValue(rec, call, sample, null, params));
    assertTrue(VcfFormatField.GL.hasValue(rec, call, sample, null, params));
    sample = VariantOutputVcfFormatterTest.createSample(Ploidy.DIPLOID, "A:G", false, 20.0, VariantSample.DeNovoStatus.IS_DE_NOVO, 10.0);
    sample.setHoeffdingAlleleBalanceHet(1.0);
    assertTrue(VcfFormatField.ABP.hasValue(rec, call, sample, null, params));
    sample = VariantOutputVcfFormatterTest.createSample(Ploidy.DIPLOID, "A:G", false, 20.0, VariantSample.DeNovoStatus.IS_DE_NOVO, 10.0);
    sample.setHoeffdingAlleleBalanceHom(1.0);
    assertTrue(VcfFormatField.ABP.hasValue(rec, call, sample, null, params));
    sample.setStatisticsString("");
    assertFalse(VcfFormatField.RS.hasValue(rec, call, sample, null, params));
  }

  public void testUpdateRecord() {
    final String[] sampleNames = {"SAMPLE"};
    final VariantParams params = VariantParams.builder().expectedCoverage(new DummyCoverageThreshold()).create();
    final double[] measures = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.1, 0.1};
    final GenotypeMeasure measure = new ArrayGenotypeMeasure(SimplePossibility.SINGLETON, measures, new HypothesesSnp(SimplePossibility.SINGLETON, GenomePriorParams.builder().create(), false, 0));
    final VariantSample sample = new VariantSample(Ploidy.DIPLOID, "A:G", false, measure, VariantSample.DeNovoStatus.IS_DE_NOVO, 10.0);
    sample.setCoverage(10);
    sample.setCoverageCorrection(5.5);
    sample.setHoeffdingAlleleBalanceHet(4.0);
    sample.setHoeffdingStrandBiasAllele1(0.5);
    sample.setHoeffdingStrandBiasAllele2(1.0);
    sample.setHoeffdingUnmatedBiasAllele1(0.3);
    sample.setHoeffdingUnmatedBiasAllele2(0.5);
    sample.setHoeffdingReadPositionBias(5.0);
    sample.setAmbiguityRatio(5.0);
    sample.setStats(new StatisticsSnp(DescriptionSnp.SINGLETON));
    sample.setStatisticsString("Q");
    sample.setPlacedUnmappedRatio(0.42);
    final Map<Set<String>, Double> like = new HashMap<>();
    like.put(Collections.singleton("A"), 0.2);
    like.put(VariantSample.pairSet("G", "A"), 0.5);
    like.put(Collections.singleton("G"), 0.2);
    sample.setGenotypeLikelihoods(like);
    Variant call = new Variant(new VariantLocus("ref", 2, 3, "A", 'G'), sample);
    VcfRecord rec = new VcfRecord("ref", 2, "A");
    rec.setNumberOfSamples(sampleNames.length);
    for (VcfFormatField field : EnumSet.range(VcfFormatField.GT, VcfFormatField.PD)) {
      field.updateRecord(rec, call, sampleNames, params, false);
    }
    assertEquals("ref\t3\t.\tA\tG\t.\t.\t.\tGT:VA:DP:DPR:RE:AR:RQ:GQ:RP:DN:DNP:ABP:SBP:RPB:PPB:AQ:PUR:RS:ADE:AD:SSC:SS:GL:GQD:ZY:PD\t0/1:.:10:14.286:5.500:5.000:10.4:1:-0.7:Y:43:4.00:1.00:5.00:0.50:0.000,0.000:0.42:Q:0.0,0.0:0,0:4.3:2:-0.53,-0.39,-0.53:0.100:e:d", rec.toString());

    sample.setHoeffdingAlleleBalanceHom(3.0);
    rec = new VcfRecord("ref", 2, "A");
    rec.setNumberOfSamples(sampleNames.length);
    for (VcfFormatField field : EnumSet.range(VcfFormatField.GT, VcfFormatField.PD)) {
      field.updateRecord(rec, call, sampleNames, params, false);
    }
    assertEquals("ref\t3\t.\tA\tG\t.\t.\t.\tGT:VA:DP:DPR:RE:AR:RQ:GQ:RP:DN:DNP:ABP:SBP:RPB:PPB:AQ:PUR:RS:ADE:AD:SSC:SS:GL:GQD:ZY:PD\t0/1:.:10:14.286:5.500:5.000:10.4:1:-0.7:Y:43:3.00:1.00:5.00:0.50:0.000,0.000:0.42:Q:0.0,0.0:0,0:4.3:2:-0.53,-0.39,-0.53:0.100:e:d", rec.toString());

    rec = new VcfRecord("ref", 2, "A");
    rec.setNumberOfSamples(sampleNames.length);
    call = new Variant(new VariantLocus("ref", 2, 3, "A", 'G'), (VariantSample) null);
    for (VcfFormatField field : EnumSet.range(VcfFormatField.GT, VcfFormatField.PD)) {
      field.updateRecord(rec, call, sampleNames, params, false);
    }

    assertEquals("ref\t3\t.\tA\t.\t.\t.\t.\tGT:VA:DP:DPR:RE:AR:RQ:GQ:RP:DN:DNP:ABP:SBP:RPB:PPB:AQ:PUR:RS:ADE:AD:SSC:SS:GL:GQD:ZY:PD\t.", rec.toString());

    call = new Variant(new VariantLocus("ref", 2, 3, "A", 'G'), sample);
    call.addFilter(VariantFilter.FAILED_COMPLEX);
    rec = new VcfRecord("ref", 2, "A");
    rec.setNumberOfSamples(sampleNames.length);
    for (VcfFormatField field : EnumSet.range(VcfFormatField.GT, VcfFormatField.PD)) {
      field.updateRecord(rec, call, sampleNames, params, false);
    }

    // This is a bit dumb the final 3 format fields haven't been populated not entirely sure why
    assertEquals("ref\t3\t.\tA\t.\t.\t.\t.\tGT:VA:DP:DPR:RE:AR:RQ:GQ:RP:DN:DNP:ABP:SBP:RPB:PPB:AQ:PUR:RS:ADE:AD:SSC:SS:GL:GQD:ZY:PD\t./.:.:10:14.286:5.500:5.000:10.4:.:.:.:43:3.00:1.00:5.00:0.50:0.000:0.42:Q:0.0:0:4.3:2:0.00", rec.toString());
  }

  public void testDenovoUpdateRecord() {
    final VariantParams params = VariantParams.builder().create();
    final VariantSample sampleA = VariantOutputVcfFormatterTest.createSample(Ploidy.DIPLOID, "A:G", false, 20.0, VariantSample.DeNovoStatus.IS_DE_NOVO, 10.0);
    final VariantSample sampleB = VariantOutputVcfFormatterTest.createSample(Ploidy.DIPLOID, "A:A", true, 20.0, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0);
    final Variant call = new Variant(new VariantLocus("ref", 2, 3, "A", 'G'), sampleA, sampleB);
    final VcfRecord rec = new VcfRecord("ref", 2, "A");
    rec.setNumberOfSamples(2);
    VcfFormatField.GT.updateRecordSample(rec, call, sampleA, null, params, false);
    VcfFormatField.DN.updateRecordSample(rec, call, sampleA, null, params, false);
    VcfFormatField.GT.updateRecordSample(rec, call, sampleB, null, params, false);
    VcfFormatField.DN.updateRecordSample(rec, call, sampleB, null, params, false);
    assertEquals("ref\t3\t.\tA\tG\t.\t.\t.\tGT:DN\t0/1:Y\t0/0:N", rec.toString());
  }

  public void testGenotypeUpdateRecord() {
    final VariantParams params = VariantParams.builder().create();
    final VariantSample[] samples = {null,
        VariantOutputVcfFormatterTest.createSample(Ploidy.DIPLOID, "A:G", false, 10.0, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0),
        VariantOutputVcfFormatterTest.createSample(Ploidy.HAPLOID, "A", true, 11.0, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0),
        VariantOutputVcfFormatterTest.createSample(Ploidy.DIPLOID, "A:A", true, 12.0, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0),
        VariantOutputVcfFormatterTest.createSample(Ploidy.HAPLOID, "G", false, 13.0, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0),
        VariantOutputVcfFormatterTest.createSample(Ploidy.DIPLOID, "G:G", false, 14.0, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0),
        VariantOutputVcfFormatterTest.createSample(Ploidy.HAPLOID, null, true, 15.0, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0),
        VariantOutputVcfFormatterTest.createSample(Ploidy.DIPLOID, null, true, 16.0, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0),
    };
    final Variant call = new Variant(new VariantLocus("ref", 2, 3, "A", 'C'), samples);
    final VcfRecord rec = new VcfRecord("ref", 2, "A");
    rec.setNumberOfSamples(samples.length);
    for (VariantSample sample : samples) {
      VcfFormatField.GT.updateRecordSample(rec, call, sample, null, params, false);
    }

    assertEquals("ref\t3\t.\tA\tG\t.\t.\t.\tGT\t.\t0/1\t0\t0/0\t1\t1/1\t.\t./.", rec.toString());
  }

  public void testAddAltCall() {
    final VcfRecord rec = new VcfRecord("ref", 2, "A");
    assertEquals(0, VcfFormatField.addAltAllele("A", "A", null, rec));
    assertEquals("ref\t3\t.\tA\t.\t.\t.\t.", rec.toString());
    assertEquals(0, VcfFormatField.addAltAllele("", "", null, rec));
    assertEquals("ref\t3\t.\tA\t.\t.\t.\t.", rec.toString());
    assertEquals(1, VcfFormatField.addAltAllele("C", "A", null, rec));
    assertEquals("ref\t3\t.\tA\tC\t.\t.\t.", rec.toString());
    assertEquals(1, VcfFormatField.addAltAllele("C", "A", null, rec));
    assertEquals("ref\t3\t.\tA\tC\t.\t.\t.", rec.toString());
    assertEquals(2, VcfFormatField.addAltAllele("C", "A", 'G', rec));
    assertEquals("ref\t3\t.\tA\tC,GC\t.\t.\t.", rec.toString());
    assertEquals(3, VcfFormatField.addAltAllele("", "A", 'G', rec));
    assertEquals("ref\t3\t.\tA\tC,GC,G\t.\t.\t.", rec.toString());
  }


  public void testQa() {
    final VcfRecord record = new VcfRecord("foo", 1, "C");
    record.addAltCall("T");
    record.addFormatAndSample("GT", "0/1");
    final Calibrator calibrator = new Calibrator(new Covariate[]{new CovariateSequence(), new CovariateReadGroup()}, new ReferenceRegions());
    final CalibratedPerSequenceExpectedCoverage expectedCoverage = new CalibratedPerSequenceExpectedCoverage(calibrator, new HashMap<>(), new HashMap<>(), new RegionRestriction("foo:1+1000")) {
      @Override
      public double expectedCoverage(String sequenceName, String sampleName) {
        return 20;
      }
    };

    final VariantParams params = new VariantParamsBuilder().expectedCoverage(expectedCoverage).create();
    final VariantSample sample = new VariantSample(Ploidy.DIPLOID, "Sample", false, new MockGenotypeMeasure(0.1), VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0);
    sample.setVariantAllele("A");
    // Override methods because it's easier than attempting to increment...
    final StatisticsSnp stats = new StatisticsSnp(DescriptionSnp.SINGLETON) {
      @Override
      public AlleleStatisticsInt counts() {
        return new AlleleStatisticsInt(DescriptionSnp.SINGLETON) {
          @Override
          public double qa(int index) {
            return 20;
          }
        };
      }
    };
    sample.setStats(stats);
    final Variant call = new Variant(new VariantLocus("foo", 1, 2), 0, sample);
    final String[] sampleNames = {"Sample"};
    record.setNumberOfSamples(sampleNames.length);
    VA.updateRecord(record, call, sampleNames, params, false);
    AQ.updateRecord(record, call, sampleNames, params, false);

    assertTrue(QA.hasValue(record, call, sample, sampleNames[0], params));

    QA.updateRecord(record, call, sampleNames, params, false);
    assertEquals("20.000", record.getFormat("QA").get(0));
  }
}
