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

package com.rtg.variant.bayes.multisample.population;

import java.io.IOException;
import java.io.PrintStream;

import com.rtg.util.InvalidParamsException;
import com.rtg.util.Utils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.Talkback;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.test.NanoRegression;

import junit.framework.TestCase;

/**
 */
public class ConvergencePlotTest extends TestCase {

  private NanoRegression mNano = null;

  @Override
  public void setUp() {
    Diagnostic.setLogStream();
    mNano = new NanoRegression(this.getClass(), false);
  }

  @Override
  public void tearDown() throws Exception {
    // clear the module name so later tests don't report SlimException to the
    // Talkback system
    Talkback.setModuleName(null);
    try {
      mNano.finish();
    } finally {
      mNano = null;
    }
  }

  private void checkVaryCoverage(final Estimator estimator, final int samples, final String id) throws InvalidParamsException, IOException {
    final MemoryPrintStream mps = new MemoryPrintStream();
    new ConvergencePlot(samples, ChrType.AUTOSOMAL, new double[] {0.33, 0.67, 0.0, 0.0}, 0.05, estimator, mps.printStream()).varyCoverage(20);
    final String str = mps.toString();
    //System.err.println(str);
    mNano.check("vc_" + id, str);
  }

  private void checkVaryCoverage(final int samples, String id) throws InvalidParamsException, IOException {
    checkVaryCoverage(new NullEstimator(), samples, "null" + id);
    checkVaryCoverage(new HwEstimator(), samples, "haploid" + id);
  }

  public void testVcLarge() throws InvalidParamsException, IOException {
    checkVaryCoverage(200, "_lge");
  }

  public void testVcMedium() throws InvalidParamsException, IOException {
    checkVaryCoverage(25, "_med");
  }

  public void testVcSmall() throws InvalidParamsException, IOException {
    checkVaryCoverage(5, "_small");
  }

  private void checkIterate(final Estimator estimator, final int samples, final String id) throws InvalidParamsException, IOException {
    final MemoryPrintStream mps = new MemoryPrintStream();
    final PrintStream ps = mps.printStream();
    final ConvergencePlot convergencePlot = new ConvergencePlot(samples, ChrType.AUTOSOMAL, new double[] {0.33, 0.67, 0.0, 0.0}, 0.05, estimator, ps);
    final int totalCalls = 200;
    final int maxCoverage = 20;
    final double average = convergencePlot.iterate(maxCoverage, totalCalls);
    final String exc = Utils.realFormat(average, 3);
    //System.err.println(exc);
    ps.println(exc);
    final String str = mps.toString();
    //System.err.println(str);
    mNano.check("it_" + id, str);
  }

  private void checkIterate(final int samples, String id) throws InvalidParamsException, IOException {
    checkIterate(new NullEstimator(), samples, "null" + id);
    checkIterate(new HwEstimator(), samples, "haploid" + id);
  }

  public void testItLarge() throws InvalidParamsException, IOException {
    checkIterate(200, "_lge");
  }

  public void testItMed() throws InvalidParamsException, IOException {
    checkIterate(25, "_med");
  }

  public void testItSmall() throws InvalidParamsException, IOException {
    checkIterate(5, "_small");
  }

}
