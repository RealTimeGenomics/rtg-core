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
package com.rtg.variant.bayes.multisample.family;


import java.io.BufferedReader;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.EnumSet;
import java.util.List;

import com.rtg.relation.Family;
import com.rtg.relation.RelationshipsFileParser;
import com.rtg.util.TestUtils;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.GenomePriorParamsBuilder;
import com.rtg.variant.StaticThreshold;
import com.rtg.variant.Variant;
import com.rtg.variant.Variant.VariantFilter;
import com.rtg.variant.VariantOutputLevel;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantParamsBuilder;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.multisample.HaploidDiploidHypotheses;
import com.rtg.variant.bayes.snp.HypothesesNone;
import com.rtg.variant.bayes.snp.HypothesesPrior;
import com.rtg.variant.format.VariantOutputVcfFormatter;
import com.rtg.variant.format.VcfInfoField;

import junit.framework.TestCase;

/**
 */
public class DiseasedFamilyCallerTest extends TestCase {

  private VariantOutputVcfFormatter makeFormatter(int numSamples) {
    final List<String> names = new ArrayList<>();
    for (int i = 0; i < numSamples; ++i) {
      names.add("g" + i);
    }
    final VariantOutputVcfFormatter formatter = new VariantOutputVcfFormatter(names.toArray(new String[0]));
    formatter.addExtraInfoFields(EnumSet.of(VcfInfoField.DPS));
    return formatter;
  }

  @SuppressWarnings("unchecked")
  public void testComparison() throws Exception {
    final GenomePriorParams params = new GenomePriorParamsBuilder().create();
    final Family f = Family.getFamily(RelationshipsFileParser.load(new BufferedReader(new StringReader(DiseasedFamilyPosteriorTest.RELATIONS))));
    assertEquals("father", f.getFather());
    assertEquals("mother", f.getMother());
    final List<ModelInterface<?>> models =  FamilyCallerTest.buildFamily(params, 1, "AAAAAAAAAAA", "AAAAAAATTTTTTTT", "AAAAAAAAAAAAAAA", "AAAAAAAATTTTTTTTTTT");
    final VariantParams vParams = new VariantParamsBuilder().maxCoverageFilter(new StaticThreshold(15)).genomePriors(params).create();
    final DiseasedFamilyCaller fc = new DiseasedFamilyCaller(vParams, f, 0.95);
    final byte[] ref = new byte[21];
    ref[19] = 3;
    ref[20] = 1;
    final Variant v = fc.makeCall("foo", 20, 21, ref, models, new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, null, (HypothesesPrior<Description>) models.get(0).hypotheses()));
    assertNotNull(v);
    final String[] outstr = makeFormatter(4).formatCall(v).split("\t");
    assertEquals(Arrays.toString(outstr), 13, outstr.length);
    final String expected = ("foo 21 . A T . OC DPS=41.0;DP=60;CT=0 GT:DP:RE:AR:GQ:ABP:SBP:RPB:AQ:PUR:RS:AD "
        + "0/0:11:2.090:0.000:64:0.00:23.89:0.00:79.337,0.000:0.00:A,11,2.090:11,0 "
        + "0/1:15:2.850:0.000:136:0.14:17.37:0.00:50.487,57.700:0.00:A,7,1.330,T,8,1.520:7,8 "
        + "0/0:15:2.850:0.000:105:0.00:32.57:0.00:108.187,0.000:0.00:A,15,2.850:15,0 "
        + "0/1:19:3.610:0.000:137:1.03:23.89:0.00:57.700,79.337:0.00:A,8,1.520,T,11,2.090:8,11\n").replaceAll(" ", "\t");

    assertEquals(expected, makeFormatter(4).formatCall(v));
    assertTrue(v.isFiltered(VariantFilter.COVERAGE));
  }

  @SuppressWarnings("unchecked")
  public void testComparisonOvercoverage() throws Exception {
    final GenomePriorParams params = new GenomePriorParamsBuilder().create();
    final Family f = Family.getFamily(RelationshipsFileParser.load(new BufferedReader(new StringReader(DiseasedFamilyPosteriorTest.RELATIONS))));
    assertEquals("father", f.getFather());
    assertEquals("mother", f.getMother());
    final List<ModelInterface<?>> models =  FamilyCallerTest.buildFamily(params, 1, "AAAAAAAAAAA", "AAAAAAATTTTTTTT", "AAAAAAAAAAAAAAA", "AAAAAAAATTTTTTTTTTT");
    final VariantParams vParams = new VariantParamsBuilder().maxCoverageFilter(new StaticThreshold(5)).callLevel(VariantOutputLevel.ALL).genomePriors(params).create();
    final DiseasedFamilyCaller fc = new DiseasedFamilyCaller(vParams, f, 0.95);
    final byte[] ref = new byte[21];
    ref[19] = 3;
    ref[20] = 1;
    final Variant v = fc.makeCall("foo", 20, 21, ref, models, new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, null, (HypothesesPrior<Description>) models.get(0).hypotheses()));
    assertNotNull(v);
    assertTrue(v.isFiltered(VariantFilter.COVERAGE));
    final String call = makeFormatter(4).formatCall(v).replace('\t', ' ');
    TestUtils.containsAll(call, "foo 21 . A T . OC DPS=41.0");
  }

  @SuppressWarnings("unchecked")
  public void testComparisonNoExplanation() throws Exception {
    final GenomePriorParams params = new GenomePriorParamsBuilder().create();
    final Family f = Family.getFamily(RelationshipsFileParser.load(new BufferedReader(new StringReader(DiseasedFamilyPosteriorTest.RELATIONS))));
    assertEquals("father", f.getFather());
    assertEquals("mother", f.getMother());
    final List<ModelInterface<?>> models =  FamilyCallerTest.buildFamily(params, 1, "AAAAAAAAAAA", "AAAAAAATTTTTTTT", "AAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAAAAA");
    final VariantParams vParams = new VariantParamsBuilder().maxCoverageFilter(new StaticThreshold(15)).genomePriors(params).create();
    final DiseasedFamilyCaller fc = new DiseasedFamilyCaller(vParams, f, 0.95);
    final byte[] ref = new byte[21];
    ref[19] = 3;
    ref[20] = 1;
    final Variant v = fc.makeCall("foo", 20, 21, ref, models, new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, null, (HypothesesPrior<Description>) models.get(0).hypotheses()));
    assertNull(v);
  }

  @SuppressWarnings("unchecked")
  public void testComparisonNoExplanationAllOutput() throws Exception {
    final GenomePriorParams params = new GenomePriorParamsBuilder().create();
    final Family f = Family.getFamily(RelationshipsFileParser.load(new BufferedReader(new StringReader(DiseasedFamilyPosteriorTest.RELATIONS))));
    assertEquals("father", f.getFather());
    assertEquals("mother", f.getMother());
    final List<ModelInterface<?>> models =  FamilyCallerTest.buildFamily(params, 1, "AAAAAAAAAAA", "AAAAAAATTTTTTTT", "AAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAAAAA");
    final VariantParams vParams = new VariantParamsBuilder()
    .callLevel(VariantOutputLevel.ALL)
    .maxCoverageFilter(new StaticThreshold(25))
    .genomePriors(params)
    .create();
    final DiseasedFamilyCaller fc = new DiseasedFamilyCaller(vParams, f, 0.95);
    final byte[] ref = new byte[21];
    ref[19] = 3;
    ref[20] = 1;
    final Variant v = fc.makeCall("foo", 20, 21, ref, models, new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, null, (HypothesesPrior<Description>) models.get(0).hypotheses()));
    assertNotNull(v);
    final String formatted = makeFormatter(4).formatCall(v);
    assertEquals(13, formatted.split("\t").length);
    assertTrue(formatted, makeFormatter(4).formatCall(v).contains("foo\t21\t.\tA\tT\t.\tPASS\tDPS=38.3;DP=59\tGT:DP:RE:AR:GQ:ABP:SBP:RPB:AQ:PUR:RS:AD\t0/0"));
  }

  @SuppressWarnings("unchecked")
  public void testShortCircuitPass() {
    final GenomePriorParams params = new GenomePriorParamsBuilder().create();
    final List<ModelInterface<?>> models =  FamilyCallerTest.buildFamily(params, 1, "AAAA", "AAAA", "AAAA", "AAAA");
    assertTrue(DiseasedFamilyCaller.shortCircuit(models, new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, null, (HypothesesPrior<Description>) models.get(0).hypotheses())));
  }

  @SuppressWarnings("unchecked")
  public void testShortCircuitNotAgree() {
    final GenomePriorParams params = new GenomePriorParamsBuilder().create();
    final List<ModelInterface<?>> models =  FamilyCallerTest.buildFamily(params, 1, "AAAA", "AAAA", "GGGG", "AAAA");
    assertFalse(DiseasedFamilyCaller.shortCircuit(models, new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, null, (HypothesesPrior<Description>) models.get(0).hypotheses())));
  }

}

