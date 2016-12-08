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
package com.rtg.simulation.reads;

import java.io.IOException;
import java.util.Arrays;

import com.rtg.reader.SdfId;
import com.rtg.util.ChiSquared;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.test.NanoRegression;
import com.rtg.variant.AbstractMachineErrorParams;

import junit.framework.TestCase;

/**
 * Test class.
 */
public abstract class AbstractMachineTest extends TestCase {

  // 0.9999 corresponds to expected failure of 1 in 10000 runs
  private static final double CONFIDENCE_LEVEL = 0.9999;

  static final int FRAGMENT_LENGTH = 1000;

  protected abstract AbstractMachineErrorParams getPriors() throws IOException, InvalidParamsException;

  protected abstract Machine getMachine(final long seed) throws IOException, InvalidParamsException;

  protected NanoRegression mNano;

  @Override
  public void setUp() {
    Diagnostic.setLogStream();
    mNano = new NanoRegression(this.getClass());
  }

  @Override
  public void tearDown() throws IOException {
    try {
      mNano.finish();
    } finally {
      mNano = null;
    }
  }


  private void checkIsProbabilityDistribution(final String name, final double[] p) {
    double s = 0;
    for (final double v : p) {
      s += v;
    }
    assertEquals(name + " " + Arrays.toString(p), 1, s, 1E-5);
  }

  void checkDiscreteDistribution(final String name, final double[] expected, final int[] actual, final double slop) {
    // Slop should be 1 for a standard test
    checkIsProbabilityDistribution(name, expected);
    assertTrue(actual.length >= expected.length);
    for (int k = expected.length; k < Math.min(2 * expected.length, actual.length); ++k) {
      // Should not be many events longer than the priors
      assertTrue(actual[k] < 10);
    }
    for (int k = 2 * expected.length; k < actual.length; ++k) {
      // Should not be no events longer than the 2 * priors length (for reasonable priors)
      assertEquals(name + " event longer than priors, length = " + k + " max-allowed=" + (expected.length - 1), 0, actual[k]);
    }
    double sum = 0;
    for (final int a : actual) {
      sum += a;
    }
    double chi = 0;
    int dofCorrection = 1;
    for (int k = 0; k < expected.length; ++k) {
      final double e = sum * expected[k];
      if (e >= 1) {
        final double g = actual[k] - e;
        final double yatesCorrection = Math.abs(g) - 0.5;
        chi += yatesCorrection * yatesCorrection / e;
        //chi += g * g / e;
      } else {
        ++dofCorrection;
        if (e == 0) {
          assertEquals(name, 0, actual[k]);
        }
      }
    }
    final int dof = Math.max(1, expected.length - dofCorrection);
    final double lower = ChiSquared.chi(dof, 1 - CONFIDENCE_LEVEL);
    final double upper = ChiSquared.chi(dof, CONFIDENCE_LEVEL);
    if (chi > slop * upper || chi * slop < lower) {
      final double[] a = new double[expected.length];
      for (int k = 0; k < a.length; ++k) {
        a[k] = actual[k] / sum;
      }
      final String ex = "Expect(" + name + "): "
      + lower + " < " + chi + " < " + upper + " events=" + sum
      + "\n" + Arrays.toString(expected) + "\n" + Arrays.toString(a);
      fail(ex);
    }
  }

  private class StatsReadWriter implements ReadWriter {

    private final int[] mMatches = new int[FRAGMENT_LENGTH];
    private final int[] mMismatches = new int[FRAGMENT_LENGTH];
    private final int[] mInserts = new int[FRAGMENT_LENGTH];
    private final int[] mDeletions = new int[FRAGMENT_LENGTH];
    private int mTotalMismatchEvents = 0;
    private int mTotalInsertEvents = 0;
    private int mTotalDeleteEvents = 0;
    private int mTotalMatchEvents = 0;
    private int mTotal = 0;
    private int mForward = 0;
    private int mLeft = 0;
    private int mRight = 0;

    private void augment(final String name, final byte[] data, final byte[] qual, final int length) {
      ++mTotal;
      if (name.contains("/F/")) {
        ++mForward;
      }
      final int lastSlash = name.lastIndexOf('/');
      assertTrue(lastSlash != -1);
      assertTrue(lastSlash != name.length() - 1);
      final String preCigar = name.substring(lastSlash + 1);
      final int colon = preCigar.indexOf(':');
      final String cigar = colon == -1 ? preCigar : preCigar.substring(0, colon);
      // Unroll cigar, incrementing appropriate events
      for (int k = 0; k < cigar.length(); ++k) {
        int n = 0;
        char c;
        while (Character.isDigit(c = cigar.charAt(k))) {
          n *= 10;
          n += c - '0';
          ++k;
        }
        switch (c) {
          case '.': // substitute for = since = not permitted in SAM names
            mMatches[n]++;
            mTotalMatchEvents += n; // for rate test, consider each unchanged position as a separate event
            break;
          case 'X':
            mMismatches[n]++;
            ++mTotalMismatchEvents;
            break;
          case 'I':
            mInserts[n]++;
            ++mTotalInsertEvents;
            break;
          case 'D':
            mDeletions[n]++;
            ++mTotalDeleteEvents;
            break;
          case 'N':
          case 'B':
            break;
          default:
            fail("Bad ciggie:" + name);
        }
      }
    }

    @Override
    public void writeRead(final String name, final byte[] data, final byte[] qual, final int length) throws IOException {
      augment(name, data, qual, length);
    }

    @Override
    public void writeLeftRead(final String name, final byte[] data, final byte[] qual, final int length) throws IOException {
      ++mLeft;
      augment(name, data, qual, length);
    }

    @Override
    public void writeRightRead(final String name, final byte[] data, final byte[] qual, final int length) throws IOException {
      ++mRight;
      augment(name, data, qual, length);
    }

    @Override
    public void identifyTemplateSet(SdfId... templateIds) { }

    @Override
    public void identifyOriginalReference(SdfId referenceId) { }

    void performStatisticalTests(final boolean paired) throws IOException, InvalidParamsException {
      // Some fundamental constraints
      assertEquals(0, mInserts[0]);
      assertEquals(0, mDeletions[0]);
      assertEquals(0, mMismatches[0]);

      if (paired) {
        assertEquals(mLeft, mRight);
        assertEquals(mTotal, mLeft + mRight);
      } else {
        assertEquals(0, mLeft);
        assertEquals(0, mRight);
      }

      // Check forward and reverse complement are equiprobable
      checkDiscreteDistribution("RC", new double[] {0.5, 0.5}, new int[] {mForward, mTotal - mForward}, 2);

      final AbstractMachineErrorParams priors = getPriors();
      final double mnpRate = priors.errorMnpEventRate();
      final double insRate = priors.errorInsEventRate();
      final double delRate = priors.errorDelEventRate();
      final double matchRate = 1 - mnpRate - insRate - delRate;
      checkDiscreteDistribution("Rates", new double[] {matchRate, mnpRate, insRate, delRate}, new int[] {mTotalMatchEvents, mTotalMismatchEvents, mTotalInsertEvents, mTotalDeleteEvents}, 5);
      checkDiscreteDistribution("MNP", priors.errorMnpDistribution(), mMismatches, getSlop());
      checkDiscreteDistribution("Del", priors.errorDelDistribution(), mDeletions, getSlop());
      checkDiscreteDistribution("Ins", priors.errorInsDistribution(), mInserts, getSlop());
    }

    @Override
    public void close() throws IOException {
    }

    @Override
    public int readsWritten() {
      return mTotal;
    }

  }

  /**
   * Slop permitted in the statistical test.  Should be 1 if the machine is
   * well modelled, but may need to be higher if the modelling is less
   * accurate.
   *
   * @return an <code>int</code> value
   */
  protected double getSlop() {
    return 1.0;
  }

  public void testDistributions() throws Exception {
    try (StatsReadWriter w = new StatsReadWriter()) {
      //      final Machine m = getMachine(System.currentTimeMillis());
      final Machine m = getMachine(42);
      m.setReadWriter(w);
      final byte[] frag = new byte[FRAGMENT_LENGTH];
      Arrays.fill(frag, (byte) 1);
      for (int k = 0; k < 10000; ++k) {
        m.processFragment("d/", 0, frag, frag.length / 2);
      }
      w.performStatisticalTests(m.isPaired());
    }
  }
}
