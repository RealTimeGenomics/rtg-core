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
package com.rtg.variant.bayes.complex;

import java.io.IOException;

import com.rtg.util.InvalidParamsException;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.test.FileHelper;
import com.rtg.variant.GenomePriorParams;

import junit.framework.TestCase;

/**
 */
public class InsertionPriorTest extends TestCase {

  public void test() throws InvalidParamsException, IOException {
    final GenomePriorParams pp = GenomePriorParams.builder().genomePriors("testhumanprior").create();
    final InsertionPrior ip = new InsertionPrior(pp);
    ip.integrity();
  }

  public void testKnownPriors() {
    final GenomePriorParams pp = GenomePriorParams.builder()
        .create();
    final InsertionPrior ip = new InsertionPrior(pp);
    ip.integrity();
  }

  public void testKnownPriors00() {
    final GenomePriorParams pp = GenomePriorParams.builder()
        .create();
    final InsertionPrior ip = new InsertionPrior(pp);
    ip.integrity();
  }

  public void testKnownPriors1() {
    final GenomePriorParams pp = GenomePriorParams.builder()
        .create();
    final InsertionPrior ip = new InsertionPrior(pp);
    ip.integrity();
  }

  public void testKnownPriors2() {
    final GenomePriorParams pp = GenomePriorParams.builder()
        .create();
    final InsertionPrior ip = new InsertionPrior(pp);
    ip.integrity();
  }

  public void testKnownPriorsn() {
    final GenomePriorParams pp = GenomePriorParams.builder()
        .create();
    final InsertionPrior ip = new InsertionPrior(pp);
    ip.integrity();
  }

  public void testKappaPriors() throws InvalidParamsException, IOException {
    Diagnostic.setLogStream();
    //see kappaimplementation.xls for details
    final GenomePriorParams pp = GenomePriorParams.builder()
        .genomeIndelDistribution(new double[] {0.6, 0.4})
        .genomeIndelEventRate(0.05)
        .genomeIndelEventFraction(0.4)
        .genomeIndelLengthDecay(0.2)
        .create();
    final InsertionPrior ip = new InsertionPrior(pp);
    ip.integrity();
    assertEquals(0.4285714, ip.heteroPriorKappa().pi(0), 1e-4);
    assertEquals(0.3571429, ip.heteroPriorKappa().piSum(1), 1e-4);
    assertEquals(0.9047619, ip.homoPriorKappa().pi(0), 1e-4);
    assertEquals(0.0595238, ip.homoPriorKappa().piSum(1), 1e-4);
    assertEquals(0.0476190, ip.homoPriorKappa().pi(1), 1e-4);
    assertEquals(6.0952381e-7, ip.homoPriorKappa().pi(8), 1e-10);

    //System.out.println(ip.toString());
    final String exp = FileHelper.resourceToString("com/rtg/variant/bayes/complex/resources/insertionpriortest0.txt");
    assertTrue(TestUtils.sameLines(exp, ip.toString(), false));
  }

  public void testDefaultHuman() {
    final GenomePriorParams pp = GenomePriorParams.builder().create();
    final InsertionPrior ip = new InsertionPrior(pp);
    ip.integrity();
  }

  public void testPriorDistribution() {
    final GenomePriorParams pp = GenomePriorParams.builder().create();
    final InsertionPrior ip = new InsertionPrior(pp);
    ip.integrity();
  }
}
