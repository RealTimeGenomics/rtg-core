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
