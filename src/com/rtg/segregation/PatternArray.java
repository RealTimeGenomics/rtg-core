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

package com.rtg.segregation;

import java.util.Arrays;

import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Encapsulate an array of patterns.
 */
//TODO pack everything into a long
public class PatternArray extends IntegralAbstract {

  /**
   * Compare two patterns and see if they can be explained by single crossover either in father or mother.
   * @param before first pattern.
   * @param after second pattern.
   * @return indicator of whether a valid cross over.
   */
  static WhichParent diff(final PatternArray before, final PatternArray after) {
    final PatternDiff fa = diff(before, after, true);
    final PatternDiff mo = diff(before, after, false);
    return new WhichParent(fa, mo);
  }

  static PatternDiff diff(final PatternArray before, final PatternArray after, final boolean select) {
    final PatternDiff a = diff(before, after, select, true);
    final PatternDiff b = diff(before, after, select, false);
    if (b.noDifference() || a.noExplantion()) {
      return b;
    }
    if (a.noDifference() || b.noExplantion()) {
      return a;
    }
    throw new RuntimeException();
  }

  private static PatternDiff diff(PatternArray before, PatternArray after, boolean select, boolean flip) {
    int diffCount = 0;
    int which = -1;
    for (int i = 0; i < after.length(); ++i) {
      final String sThis = after.index(i).sString(select);
      final String sThat = before.index(i).sString(select);
      if (!sThis.equals(sThat) ^ flip) {
        ++diffCount;
        which = i;
      }
    }
    if (diffCount == 0) {
      return new PatternDiff(true, false, -1, flip);
    } else if (diffCount == 1) {
      return new PatternDiff(false, false, which, flip);
    } else {
      return new PatternDiff(false, true, -1, false);
    }
  }

  /**
   * @param before the first block being checked
   * @param after  the second block being checked
   * @param xLike  true if treating like an X chromosome
   * @return a cross over if the difference between the two blocks can be explained by a cross-over (otherwise null).
   */
  static CrossOver crossover(final PatternArray before, final PatternArray after, boolean xLike) {
    assert after.length() == before.length();
    if (!before.isUniqueAndValid(xLike) || !after.isUniqueAndValid(xLike)) {
      return null;
    }
    final WhichParent diff = diff(before, after);
    if (!diff.isValid()) {
      return null;
    }
    final CrossOver xo = PatternArray.crossoverInternal(before, after, diff);
    assert xo != null;
    return xo;
  }

  private static CrossOver crossoverInternal(final PatternArray before, PatternArray after, final WhichParent diff) {
    final boolean select = diff.isFather();
    final PatternArray consistent = after.flip(Pattern.faMo2Flip(diff.father().flipped(), diff.mother().flipped()));

    assert after.length() == before.length();
    final int beforePhaseCount = before.phaseCount(select);
    final int afterPhaseCount = consistent.phaseCount(select);
    return new CrossOver(before.length(), select, diff.child(), beforePhaseCount, afterPhaseCount, consistent);
  }

  private final Pattern[] mPatterns;

  private long[] mSignature = null;

  private int mHash;

  PatternArray(final String fa, final String mo) {
    assert fa.length() == mo.length();
    //assert fa.matches("^[01]+$") && mo.matches("^[01]+$");
    final Pattern[] patterns = new Pattern[fa.length()];
    for (int i = 0; i < fa.length(); ++i) {
      patterns[i] = Pattern.stringToPattern(fa.charAt(i), mo.charAt(i));
    }
    mPatterns = patterns;
    assert integrity();
  }

  PatternArray(final Pattern... patterns) {
    mPatterns = patterns;
    assert integrity();
  }

  //For testing
  PatternArray(final int... q) {
    final Pattern[] p = new Pattern[q.length];
    for (int i = 0; i < q.length; ++i) {
      p[i] = new Pattern(q[i]);
    }
    mPatterns = p;
    assert integrity();
  }

  Pattern index(final int i) {
    return mPatterns[i];
  }

  /**
   * @param selectFa true if to return father count, otherwise mother count.
   * @return the count of number of children in one of the phasings.
   */
  Integer phaseCount(final boolean selectFa) {
    int cnt = 0;
    for (int i = 0; i < length(); ++i) {
      final Pattern pat = index(i);
      final String str = selectFa ? pat.faString() : pat.moString();
      switch (str) {
        case "0":
          ++cnt;
          break;
        case "1":
          //do nothing
          break;
        case "?":
          return null;
        default:
          throw new RuntimeException();
      }
    }
    return cnt;
  }

  /**
   * @return number of children in one of the phasings for father.
   */
  int faCount() {
    int cnt = 0;
    for (Pattern mPattern : mPatterns) {
      if ("0".equals(mPattern.faString())) {
        ++cnt;
      }
    }
    return cnt;
  }

  /**
   * @return number of children in one of the phasings for mother.
   */
  int moCount() {
    int cnt = 0;
    for (Pattern mPattern : mPatterns) {
      if ("0".equals(mPattern.moString())) {
        ++cnt;
      }
    }
    return cnt;
  }

  /**
   * @return human readable string with phasings for father.
   */
  String faString() {
    final StringBuilder sb = new StringBuilder();
    for (Pattern mPattern : mPatterns) {
      sb.append(mPattern.faString());
    }
    return sb.toString();
  }

  /**
   * @return human readable string with phasings for mother.
   */
  String moString() {
    final StringBuilder sb = new StringBuilder();
    for (Pattern mPattern : mPatterns) {
      sb.append(mPattern.moString());
    }
    return sb.toString();
  }

  /**
   * @param xChr true if treating as X chromosome
   * @return true iff all patterns have a unique valid phasing (there will be no "?" in the output).
   */
  boolean isUniqueAndValid(boolean xChr) {
    for (Pattern mPattern : mPatterns) {
      if (!(xChr ? mPattern.isUniqueAndValidXChromosome() : mPattern.isUniqueAndValid())) {
        return false;
      }
    }
    return true;
  }

  /**
   * @return the length of the patterns (the number of children in the family).
   */
  int length() {
    return mPatterns.length;
  }

  /**
   * @param pa  pattern array to be checked.
   * @return true iff the two pattern arrays are compatible.
   */
  boolean compatible(final PatternArray pa) {
    assert mPatterns.length == pa.length();
    for (int i = 0; i < Pattern.NUMBER_FLIPS; ++i) {
      if (compatible(pa, i)) {
        return true;
      }
    }
    return false;
  }

  private boolean compatible(final PatternArray pa, final int flip) {
    assert mPatterns.length == pa.length();
    final int length = Math.min(length(), pa.length());
    for (int i = 0; i < length; ++i) {
      if (Pattern.flipIntersect(mPatterns[i], pa.mPatterns[i], flip) == null) {
        return false;
      }
    }
    return true;
  }

  /**
   * @param pa pattern array to be flipped and intersected.
   * @param flip identifies which flip to do.
   * @return the flipped and intersected pattern array - null if no valid intersection.
   */
  PatternArray flipIntersect(final PatternArray pa, final int flip) {
    assert length() == pa.length();
    final Pattern[] newPat = new Pattern[length()];
    for (int i = 0; i < length(); ++i) {
      final Pattern flp = Pattern.flipIntersect(mPatterns[i], pa.mPatterns[i], flip);
      if (flp == null) {
        return null;
      }
      newPat[i] = flp;
    }
    return new PatternArray(newPat);
  }

  //Used as a utility in testing only
  PatternArray flip(final int flip) {
    final int[] newPat = new int[length()];
    for (int i = 0; i < length(); ++i) {
      newPat[i] = mPatterns[i].flip(flip);
    }
    return new PatternArray(newPat);
  }

  @Override
  public int hashCode() {
    checkSignature();
    return mHash;
  }

  /**
   * Pack the pattern array into a long after a pattern flip.
   * @param flip the flip to use.
   * @return the packed long.
   */
  private long signature(final int flip) {
    assert Pattern.NUMBER_BITS * length() <= Long.SIZE;
    long t = 0;
    for (int i = 0; i < length(); ++i) {
      t = (t << Pattern.NUMBER_BITS) | mPatterns[i].flip(flip);
    }
    return t;
  }

  @Override
  public boolean equals(Object obj) {
    if (this == obj) {
      return true;
    }
    if (obj == null) {
      return false;
    }
    final PatternArray that = (PatternArray) obj;
    this.checkSignature();
    that.checkSignature();
    return Arrays.equals(this.mSignature, that.mSignature);
  }

  /**
   * @param obj pattern array to be tested for strict equality.
   * @return true iff this pattern has an identical pattern to obj. No allowance is made for the ambiguity and flipping of the phasing.
   */
  boolean strictEquals(Object obj) {
    if (this == obj) {
      return true;
    }
    if (obj == null) {
      return false;
    }
    final PatternArray that = (PatternArray) obj;
    return Arrays.equals(this.mPatterns, that.mPatterns);
  }

  /**
   * Lazily compute a signature that is unique given possible flips of patterns, also compute the hash from this.
   */
  private void checkSignature() {
    if (mSignature != null) {
      return;
    }
    mSignature = new long[Pattern.NUMBER_FLIPS];
    for (int i = 0; i < Pattern.NUMBER_FLIPS; ++i) {
      mSignature[i] = signature(i);
    }
    Arrays.sort(mSignature);
    final int prime = 31;
    long t = prime;
    for (final long f : mSignature) {
      t = t * prime + f;
    }
    final long tt = t * 3221225461L; //large prime < 2^32
    final long tx = tt + (tt >>> 32);
    mHash = (int) tx;
  }

  @Override
  public void toString(StringBuilder sb) {
    sb.append(" fa: ").append(faString());
    sb.append(" mo: ").append(moString());
  }

  @Override
  public boolean integrity() {
    Exam.assertNotNull(mPatterns);
    if (mSignature != null) {
      Exam.assertEquals(Pattern.NUMBER_FLIPS, mSignature.length);
      final int n = length() * Pattern.NUMBER_BITS;
      final long mask = -1L << n;
      for (int i = 0; i < Pattern.NUMBER_FLIPS; ++i) {
        Exam.assertEquals(0, mSignature[i] & mask);
      }
    }
    return true;
  }

}
