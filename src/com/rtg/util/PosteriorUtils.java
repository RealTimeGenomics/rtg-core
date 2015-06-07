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
package com.rtg.util;


/**
 * Container for posterior conversion methods. All our posteriors are represented as natural log scaled odds.
 *
 */
public final class PosteriorUtils {

  private static final double TEN_LOG_TEN = 10.0 / Math.log(10.0);

  private static final double INV_TEN_LOG_TEN = 1.0 / TEN_LOG_TEN;

  /** For values greater than this 1 + e^x == e^x for doubles. */
  private static final double MAX_EXP = 36.9;

  /** For values greater than this 10^(x / 10) == 10^(x / 10) - 1 for doubles. */
  private static final double MAX_PWR = 161.4;

  private PosteriorUtils() { }

  /**
   * @param posterior the posterior score (natural log scaled odds)
   * @return the posterior score in phred scale
   */
  public static double phredIfy(double posterior) {
    if (posterior > MAX_EXP) {
      return posterior * TEN_LOG_TEN;
    }
    return TEN_LOG_TEN * Math.log(1 + Math.exp(posterior));
  }

  /**
   * @param phredPosterior the posterior score in phred scale
   * @return the posterior score (natural log scaled odds)
   */
  public static double unphredIfy(double phredPosterior) {
    //from typing <code>y= -10*log10(1/(1+e^x)) solve for x</code> into wolfram alpha
    if (phredPosterior > MAX_PWR) {
      return INV_TEN_LOG_TEN * phredPosterior;
    }
    return Math.log(Math.pow(10, phredPosterior / 10.0) - 1);
  }

  /**
   * @param nonIdentityPosterior the non-identity posterior score (natural log scaled)
   * @return the non-identity posterior score in phred scale
   */
  public static double nonIdentityPhredIfy(double nonIdentityPosterior) {
    return phredIfy(-nonIdentityPosterior);
  }

  /**
   * @param identityPhred the phred scale non-identity posterior score
   * @return the non-identity posterior score (natural log scaled)
   */
  public static double nonIdentityUnPhredIfy(double identityPhred) {
    return -unphredIfy(identityPhred);
  }

  //  public static void main(final String[] args) {
  //    for (int i = 1580; i < 1620; i++) {
  //      final double x = i / 10.0;
  //      final double y = Math.pow(10, x / 10.0);
  //      final boolean same = y + 1 == y;
  //      System.err.println(x + " " + (same ? "-" : "+"));
  //    }
  //  }
}

