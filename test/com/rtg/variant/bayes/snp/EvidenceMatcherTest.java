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
package com.rtg.variant.bayes.snp;

import static com.rtg.util.StringUtils.LS;

import java.io.IOException;

import com.rtg.util.InvalidParamsException;
import com.rtg.util.filelinechecks.SimpleLineCheck;
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
      super(GenomePriorParams.builder().create(), true);
    }

    @Override
    protected ModelInterface<Description> makeModel(Hypotheses<Description> hyp) {
      return new Model<Description>(hyp, new StatisticsSnp(hyp.description())) {
        @Override
        public String toString() {
          return "ref=" + hypotheses().reference() + " " + statistics().toString();
        }
      };
    }

  }

  public void test1() throws InvalidParamsException, IOException {
    final ModelFactory<Description, ?> factory = new MockModelFactory();
    final byte[] template = {1, 2, 3, 4, 0};
    final EvidenceMatcher<ModelInterface<Description>> bm = new EvidenceMatcher<>(new ReferenceBasedBuffer<>(1, factory, template, 0), new EvidenceQFactory());
    assertEquals(INIT, bm.toString());
    final int phred = 63;
    final double q = VariantUtils.phredToProb(phred);

    final EvidenceQ evid = new EvidenceQ(DESC, 1, 0, 0, q, q, true, false, false, false);
    bm.match(0, evid);
    bm.match(1, evid);
    bm.match(2, evid);
    bm.match(3, evid);
    bm.match(4, evid);
    bm.match(0, evid);
    bm.match(0, evid);

    //the other calls
    bm.match(1, 0, 0, 0, phred, phred, null, bm.getStateIndex(true, false, false));
    bm.match(1, 0, 0, 1, phred, phred, null, bm.getStateIndex(true, false, false));

    assertEquals(EXP1, bm.toString());
    final VariantParams params = new VariantParamsBuilder().callLevel(VariantOutputLevel.ALL).create();
    checkCall(template, bm, params, 0, new Object[] {"G1", "1", ".", "A", "C", "154.5", "PASS", ".", "GT:DP:RE:AR:GQ:ABP:SBP:RPB:PUR:RS:AD:GL", "1:3:0.000:0.000:155:0.00:6.51:0.00:0.00:C,3,0.000:0,3:-15.45,0.00"});

    try {
      //do again at the same position
      makeCall(bm, 0, template, params);
      fail();
    } catch (final RuntimeException e) {
      assertEquals("ref=0 != base=1", e.getMessage());
    }
    assertNull(makeCall(bm, 1, template, params));
    checkCall(template, bm, params, 2, new Object[] {"G1", "3", ".", "G", "C", "25.5", "PASS", ".", "GT:DP:RE:AR:GQ:ABP:SBP:RPB:PUR:RS:AD:GL", "1:1:0.000:0.000:26:0.00:2.17:0.00:0.00:C,1,0.000:0,1:-2.55,0.00"});
    checkCall(template, bm, params, 3, new Object[] {"G1", "4", ".", "T", "C", "30.1", "PASS", ".", "GT:DP:RE:AR:GQ:ABP:SBP:RPB:PUR:RS:AD:GL", "1:1:0.000:0.000:30:0.00:2.17:0.00:0.00:C,1,0.000:0,1:-3.01,0.00"});
  }

  private Variant makeCall(EvidenceMatcher<ModelInterface<Description>> bm, int position, byte[] template, VariantParams params) {
    final ModelInterface<Description> model = bm.step(position);
    if (model == null) {
      return null;
    }
    Variant call = ModelTest.makeCalls(model, "G1", position, position + 1, template, params);
    return (call != null && call.isInteresting()) ? call : null;
  }

  private void checkCall(byte[] template, EvidenceMatcher<ModelInterface<Description>> bm, VariantParams params, int position, Object[] checks) {
    final Variant call = makeCall(bm, position, template, params);
    assertNotNull(call);
    final String bostr = FORMATTER.formatCall(call);
    //System.err.println(bostr);
    SimpleLineCheck.TAB_CHECK.lineCheck(bostr.trim(), checks);
  }
}
