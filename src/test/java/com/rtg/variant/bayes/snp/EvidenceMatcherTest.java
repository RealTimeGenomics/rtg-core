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
package com.rtg.variant.bayes.snp;

import static com.rtg.util.StringUtils.LS;

import com.rtg.util.InvalidParamsException;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.Variant;
import com.rtg.variant.VariantOutputLevel;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantParamsBuilder;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.Model;
import com.rtg.variant.bayes.ModelFactory;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.ModelTest;
import com.rtg.variant.bayes.NoAlleleBalance;
import com.rtg.variant.format.VariantOutputVcfFormatter;
import com.rtg.variant.util.VariantUtils;

import junit.framework.TestCase;

/**
 */
public class EvidenceMatcherTest extends TestCase {

  private static final VariantOutputVcfFormatter FORMATTER = new VariantOutputVcfFormatter();

  private static final String INIT = ""
      + "Matcher:" + LS
      + "Buffer length=1 base=0 current=0" + LS
      + ">>>>" + LS
      + "[0]     null" + LS
      + "" + LS
      ;
  private static final String EXP1 = ""
      + "Matcher:" + LS
      + "Buffer length=6 base=0 current=0" + LS
      + ">>>>" + LS
      + "[0]     ref=0 Counts " + LS
      + "coverage=3 correction=0.000" + LS
      + " [0]  0  0.000 [1]  3  0.000 [2]  0  0.000 [3]  0  0.000" + LS
      + "[1]     ref=1 Counts " + LS
      + "coverage=2 correction=0.000" + LS
      + " [0]  1  0.000 [1]  1  0.000 [2]  0  0.000 [3]  0  0.000" + LS
      + "[2]     ref=2 Counts " + LS
      + "coverage=1 correction=0.000" + LS
      + " [0]  0  0.000 [1]  1  0.000 [2]  0  0.000 [3]  0  0.000" + LS
      + "[3]     ref=3 Counts " + LS
      + "coverage=1 correction=0.000" + LS
      + " [0]  0  0.000 [1]  1  0.000 [2]  0  0.000 [3]  0  0.000" + LS
      + "[4]     ref=-1 Counts " + LS
      + "coverage=1 correction=0.000" + LS
      + " [0]  0  0.000 [1]  1  0.000 [2]  0  0.000 [3]  0  0.000" + LS
      + "[5]     null" + LS
      + "" + LS
      ;

  private static final DescriptionCommon DESC = DescriptionSnp.SINGLETON;
  private static class MockModelFactory extends ModelSnpFactory {

    MockModelFactory() {
      super(GenomePriorParams.builder().create(), true, new NoAlleleBalance());
    }

    @Override
    protected ModelInterface<Description> makeModel(Hypotheses<Description> hyp) {
      return new Model<Description>(hyp, new StatisticsSnp(hyp.description()), new NoAlleleBalance()) {
        @Override
        public String toString() {
          return "ref=" + hypotheses().reference() + " " + statistics().toString();
        }
      };
    }

  }

  public void test1() throws InvalidParamsException {
    final ModelFactory<Description, ?> factory = new MockModelFactory();
    final byte[] template = {1, 2, 3, 4, 0};
    final EvidenceMatcher<ModelInterface<Description>> bm = new EvidenceMatcher<>(new ReferenceBasedBuffer<>(1, factory, template, 0), new EvidenceQFactory());
    assertEquals(INIT, bm.toString());
    final int phred = 63;
    final double q = VariantUtils.phredToProb(phred);

    final EvidenceQ evid = new EvidenceQ(DESC, 1, 0, 0, q, q, true, false, true, false, false);
    bm.match(0, evid);
    bm.match(1, evid);
    bm.match(2, evid);
    bm.match(3, evid);
    bm.match(4, evid);
    bm.match(0, evid);
    bm.match(0, evid);

    //the other calls
    bm.match(1, 0, 0, 0, phred, phred, bm.getStateIndex(true, false, true, false));
    bm.match(1, 0, 0, 1, phred, phred, bm.getStateIndex(true, false, true, false));

    assertEquals(EXP1, bm.toString());
    final VariantParams params = new VariantParamsBuilder().callLevel(VariantOutputLevel.ALL).create();
    checkCall(template, bm, params, 0, "G1\t1\t.\tA\tC\t154.5\tPASS\tDP=3\tGT:DP:RE:AR:GQ:ABP:SBP:RPB:AQ:PUR:RS:AD:GL\t1:3:0.000:0.000:155:0.00:6.51:0.00:0.000,179.969:0.00:C,3,0.000:0,3:-15.45,0.00");

    try {
      //do again at the same position
      makeCall(bm, 0, template, params);
      fail();
    } catch (final RuntimeException e) {
      assertEquals("ref=0 != base=1", e.getMessage());
    }
    assertNull(makeCall(bm, 1, template, params));
    checkCall(template, bm, params, 2, "G1\t3\t.\tG\tC\t25.5\tPASS\tDP=1\tGT:DP:RE:AR:GQ:ABP:SBP:RPB:AQ:PUR:RS:AD:GL\t1:1:0.000:0.000:26:0.00:2.17:0.00:0.000,59.990:0.00:C,1,0.000:0,1:-2.55,0.00");
    checkCall(template, bm, params, 3, "G1\t4\t.\tT\tC\t30.1\tPASS\tDP=1\tGT:DP:RE:AR:GQ:ABP:SBP:RPB:AQ:PUR:RS:AD:GL\t1:1:0.000:0.000:30:0.00:2.17:0.00:0.000,59.990:0.00:C,1,0.000:0,1:-3.01,0.00");
  }

  private Variant makeCall(EvidenceMatcher<ModelInterface<Description>> bm, int position, byte[] template, VariantParams params) {
    final ModelInterface<Description> model = bm.step(position);
    if (model == null) {
      return null;
    }
    model.freeze();
    final Variant call = ModelTest.makeCalls(model, "G1", position, position + 1, template, params);
    return (call != null && call.isInteresting()) ? call : null;
  }

  private void checkCall(byte[] template, EvidenceMatcher<ModelInterface<Description>> bm, VariantParams params, int position, String checks) {
    final Variant call = makeCall(bm, position, template, params);
    assertNotNull(call);
    final String bostr = FORMATTER.formatCall(call);
    //System.err.println(bostr);
    assertEquals(checks, bostr.trim());
  }
}
