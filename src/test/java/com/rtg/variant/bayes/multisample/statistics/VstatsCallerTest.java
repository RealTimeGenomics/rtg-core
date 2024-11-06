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

package com.rtg.variant.bayes.multisample.statistics;


import static com.rtg.util.StringUtils.LS;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import com.rtg.launcher.OutputParams;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.VariantParams;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.MockHypotheses;
import com.rtg.variant.bayes.MockModel;
import com.rtg.variant.bayes.Model;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.NoAlleleBalance;
import com.rtg.variant.bayes.multisample.HaploidDiploidHypotheses;
import com.rtg.variant.bayes.snp.DescriptionSnp;
import com.rtg.variant.bayes.snp.EvidenceQ;
import com.rtg.variant.bayes.snp.HypothesesNone;
import com.rtg.variant.bayes.snp.StatisticsSnp;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 */
public class VstatsCallerTest extends TestCase {

  private static final double[] PRIORS = new double[]{0.1, 0.4, 0.35, 0.15};

  @Override
  protected void setUp() {
    Diagnostic.setLogStream();
  }

  public void test() throws Exception {
    Diagnostic.setLogStream();
    final double[] priors = PRIORS;
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final MockHypotheses<Description> hypotheses = new MockHypotheses<>(DescriptionSnp.SINGLETON, arith, true, priors, 0);
    final List<ModelInterface<?>> list = new ArrayList<>();
    final MockModel<Description> model = new  MockModel<>(hypotheses, new StatisticsSnp(hypotheses.description()), null);
    list.add(model);
    final File tempDir = FileUtils.createTempDir("vstats", "test");
    try {
      final VariantParams p = getOutputParams(tempDir);
      final VstatsCaller sc = new VstatsCaller(p);
      assertNull(sc.makeCall("test", 3, 4, new byte[] {0, 1, 2, 3, 4, 5}, list, new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, null, null)));
      sc.close();
      final File file = new File(tempDir, "haploid_statistics");
      assertTrue(file.exists());
      final String exp = ""
          + "#Coverage R X Count Total" + LS
          + "0 0 0 1 1" + LS
          ;
      assertEquals(exp, FileHelper.fileToString(file));
      final File diploid = new File(tempDir, "diploid_statistics");
      assertTrue(diploid.exists());
      final String expDiploid = ""
          + "#Coverage R X Count Total" + LS
          ;
      assertEquals(expDiploid, FileHelper.fileToString(diploid));
    } finally {
      assertTrue(FileHelper.deleteAll(tempDir));
    }
  }

  public void testNtRange() throws Exception {
    final double[] priors = PRIORS;
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final MockHypotheses<Description> hypotheses = new MockHypotheses<>(DescriptionSnp.SINGLETON, arith, true, priors, 0);
    final List<ModelInterface<?>> list = new ArrayList<>();
    final MockModel<Description> model = new  MockModel<>(hypotheses, new StatisticsSnp(hypotheses.description()), null);
    list.add(model);
    final File tempDir = FileUtils.createTempDir("vstats", "test");
    try {
      final VariantParams p = getOutputParams(tempDir);
      for (int i = 0; i < 10; ++i) {
        list.get(0).increment(new EvidenceQ(hypotheses.description(), 3, 0, 0, 0.01, 0.1, true, false, true, false, false));
      }
      final VstatsCaller sc = new VstatsCaller(p);
      for (int i = 0; i < 5; ++i) {
        assertNull(sc.makeCall("test", i, i + 1, new byte[] {0, 1, 2, 3, 4, 5}, list, new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, null, null)));
      }
      sc.close();
      final File file = new File(tempDir, "haploid_statistics");
      assertTrue(file.exists());
      final String exp = ""
          + "#Coverage R X Count Total" + LS
          + "10 0 10 3 4" + LS
          + "10 10 0 1 4" + LS
          ;
      assertEquals(exp, FileHelper.fileToString(file));
    } finally {
      assertTrue(FileHelper.deleteAll(tempDir));
    }
  }

  public void testNullModel() throws Exception {
    final File tempDir = FileUtils.createTempDir("vstats", "test");
    try {
      final VariantParams p = getOutputParams(tempDir);
      final VstatsCaller sc = new VstatsCaller(p);
      final List<ModelInterface<?>> list = new ArrayList<>();
      list.add(null);
      assertNull(sc.makeCall("test", 3, 4, new byte[]{0, 1, 2, 3, 4, 5}, list, new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, null, null)));
      sc.close();
      final File file = new File(tempDir, "haploid_statistics");
      assertTrue(file.exists());
      final String exp = ""
          + "#Coverage R X Count Total" + LS
          ;
      assertEquals(exp, FileHelper.fileToString(file));
    } finally {
      assertTrue(FileHelper.deleteAll(tempDir));
    }
  }

  public void testWithIncrements() throws Exception {
    final double[] priors = {0.1, 0.4, 0.35, 0.15};
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final MockHypotheses<Description> hypotheses = new MockHypotheses<>(DescriptionSnp.SINGLETON, arith, true, priors, 0);
    final List<ModelInterface<?>> listA = new ArrayList<>();
    listA.add(new  Model<>(hypotheses, new StatisticsSnp(hypotheses.description()), new NoAlleleBalance()));
    final List<ModelInterface<?>> listB = new ArrayList<>();
    listB.add(new  Model<>(hypotheses, new StatisticsSnp(hypotheses.description()), new NoAlleleBalance()));
    final File tempDir = FileUtils.createTempDir("vstats", "test");
    try {
      final VariantParams p = getOutputParams(tempDir);
      updateModels(hypotheses, listA, listB);
      final VstatsCaller sc = new VstatsCaller(p);
      assertNull(sc.makeCall("test", 0, 1, new byte[] {1, 1, 2, 3, 4, 5}, listA, new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, hypotheses, null)));
      assertNull(sc.makeCall("test", 3, 4, new byte[] {1, 1, 2, 3, 4, 5}, listB, new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, hypotheses, null)));
      sc.close();
      final File file = new File(tempDir, "haploid_statistics");
      assertTrue(file.exists());
      assertEquals(EXPECTED, FileHelper.fileToString(file));
    } finally {
      assertTrue(FileHelper.deleteAll(tempDir));
    }
  }

  private void updateModels(final MockHypotheses<Description> hypotheses, final List<ModelInterface<?>> listA, final List<ModelInterface<?>> listB) {
    for (int i = 0; i < 10; ++i) {
      listA.get(0).increment(new EvidenceQ(hypotheses.description(), 3, 0, 0, 0.01, 0.1, true, false, true, false, false));
    }
    for (int i = 0; i < 5; ++i) {
      listB.get(0).increment(new EvidenceQ(hypotheses.description(), 2, 0, 0, 0.01, 0.1, true, false, true, false, false));
    }
  }

  private static final String EXPECTED_EMPTY = ""
      + "#Coverage R X Count Total" + LS
      ;

  private static final String EXPECTED = ""
      + "#Coverage R X Count Total" + LS
      + "5 5 0 1 1" + LS
      + "10 0 10 1 1" + LS
      ;

  public void testWithIncrementsDiploid() throws Exception {
    final double[] priors = {0.1, 0.4, 0.35, 0.15, 0.1, 0.4, 0.35, 0.15, 0.35, 0.15};
    final MockHypotheses<Description> hypotheses = new MockHypotheses<>(DescriptionSnp.SINGLETON, SimplePossibility.SINGLETON, false, priors, 0);
    final List<ModelInterface<?>> listA = new ArrayList<>();
    listA.add(new  Model<>(hypotheses, new StatisticsSnp(hypotheses.description()), new NoAlleleBalance()));
    final List<ModelInterface<?>> listB = new ArrayList<>();
    listB.add(new  Model<>(hypotheses, new StatisticsSnp(hypotheses.description()), new NoAlleleBalance()));
    final File tempDir = FileUtils.createTempDir("vstats", "test");
    try {
      final VariantParams p = getOutputParams(tempDir);
      updateModels(hypotheses, listA, listB);
      final VstatsCaller sc = new VstatsCaller(p);
      assertNull(sc.makeCall("test", 0, 1, new byte[] {1, 1, 2, 3, 4, 5}, listA, new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, null, hypotheses)));
      assertNull(sc.makeCall("test", 3, 4, new byte[] {1, 1, 2, 3, 4, 5}, listB, new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, null, hypotheses)));
      sc.close();
      final File file = new File(tempDir, "diploid_statistics");
      assertTrue(file.exists());

      assertEquals(EXPECTED, FileHelper.fileToString(file));
      final File haploid = new File(tempDir, "haploid_statistics");
      assertTrue(haploid.exists());

      assertEquals(EXPECTED_EMPTY, FileHelper.fileToString(haploid));
    } finally {
      assertTrue(FileHelper.deleteAll(tempDir));
    }
  }

  private VariantParams getOutputParams(final File tempDir) {
    final OutputParams outputParams = new OutputParams(tempDir, false);
    return VariantParams.builder()
        .machineErrorName("illumina")
        .genomePriors(GenomePriorParams.builder().create())
        .noComplexCalls(true)
        .outputParams(outputParams)
        .create();
  }
}
