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
package com.rtg.variant.coverage;

import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicIntegerArray;
import java.util.concurrent.atomic.AtomicLongArray;

import com.rtg.util.MathUtils;

/**
 * Maintain thread-safe coverage related information for a sequence.
 */
class CoverageState {

  /** The actual bytes of the current template; read only. */
  private byte[] mTemplate;

  /** Number of reads hitting each template position, 32:32 fixed-point arithmetic. */
  private final AtomicLongArray mCoverage;

  /** Number of <code>IH=1</code> reads hitting each template position. */
  private final AtomicIntegerArray mIH1;

  /** Number of <code>IH&gt;1</code> reads hitting each template position. */
  private final AtomicIntegerArray mIHgt1;

  private final double[] mScore;

  /** True if the current template has any mappings. */
  private final AtomicBoolean mTemplateSeen = new AtomicBoolean();

  private final String mTemplateName;

  CoverageState(final String templateName, final byte[] templateBytes, final boolean additional) {
    mTemplateName = templateName;
    mTemplate = templateBytes;
    mCoverage = new AtomicLongArray(mTemplate.length);
    if (additional) {
      mIH1 = new AtomicIntegerArray(mTemplate.length);
      mIHgt1 = new AtomicIntegerArray(mTemplate.length);
      mScore = new double[mTemplate.length];
    } else {
      mIH1 = null;
      mIHgt1 = null;
      mScore = null;
    }
  }

  // Scale chosen so that IH up to 11 can be exactly represented, but
  // still small enough that only ~15 bits are used for the fraction.
  // Testing shows this is actually more accurate than IEEE double
  // arithmetic for the situation of interest
  private static final double SCALE = 16.0 * 9.0 * 5 * 7 * 11;
  private static final double INV_SCALE = 1.0 / SCALE;

  double getCoverage(final int pos) {
    final long c = mCoverage.get(pos);
    return c * INV_SCALE;
  }

  void incrementCoverage(final int pos, final double coverage) {
    final long c = MathUtils.round(coverage * SCALE);
    mCoverage.addAndGet(pos, c);
  }

  byte getTemplate(final int pos) {
    return mTemplate[pos];
  }

  int getTemplateLength() {
    return mTemplate.length;
  }

  void incrementIH(final int pos, final int ih) {
    if (ih == 1) {
      mIH1.incrementAndGet(pos);
    } else {
      mIHgt1.incrementAndGet(pos);
    }
    addToScore(pos, 1.0 / ih);
  }

  synchronized void addToScore(final int pos, final double score) {
    mScore[pos] += score;
  }

  int getIH1(final int pos) {
    return mIH1.get(pos);
  }

  int getIHgt1(final int pos) {
    return mIHgt1.get(pos);
  }

  double getScore(final int pos) {
    return mScore[pos];
  }

  void setSeen() {
    mTemplateSeen.set(true);
  }

  boolean isSeen() {
    return mTemplateSeen.get();
  }

  boolean isAdditional() {
    return mIH1 != null;
  }

  String getName() {
    return mTemplateName;
  }
}
