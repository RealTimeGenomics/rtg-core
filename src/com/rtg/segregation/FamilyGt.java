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

import java.util.Map;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.reference.Ploidy;
import com.rtg.reference.ReferenceSequence;
import com.rtg.reference.Sex;
import com.rtg.util.Pair;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * The genotypes for all members of family (parents and all children).
 */
@TestClass({"com.rtg.segregation.SegregationVcfSearchTest", "com.rtg.segregation.FamilyGtTest"})
public class FamilyGt extends IntegralAbstract {

  protected static PatternArray childPatterns(final GType fa, final GType mo, final GType[] ch) {
    final Pattern[] chp = new Pattern[ch.length];
    for (int i = 0; i < ch.length; ++i) {
      chp[i] = new Pattern(fa, mo, ch[i]);
    }
    return new PatternArray(chp);
  }

  static FamilyGt familyPloidy(final String seq, final int posn, final String[] strs, Sex[] sexes, Map<Pair<Sex, String>, ReferenceSequence> ploidyMap) throws MismatchingPloidyException {
    final GType fa = new GType(strs[0], ploidyMap.get(new Pair<>(sexes[0], seq)).effectivePloidy(posn));
    final GType mo = new GType(strs[1], ploidyMap.get(new Pair<>(sexes[1], seq)).effectivePloidy(posn));
    final GType[] ch = new GType[strs.length - 2];
    for (int i = 2; i < strs.length; ++i) {
      ch[i - 2] = new GType(strs[i], ploidyMap.get(new Pair<>(sexes[i], seq)).effectivePloidy(posn));
    }
    final boolean isXlike = fa.ploidy() == Ploidy.HAPLOID && mo.ploidy() == Ploidy.DIPLOID;
    return new FamilyGt(seq, posn, fa, mo, ch, isXlike);
  }


  private final String mSeq;
  private final int mPosn;
  private final GType mFa;
  private final GType mMo;
  private final GType[] mCh;
  private final PatternArray mChP;
  private final boolean mIsXlike;

  FamilyGt(String seq, int posn, GType fa, GType mo, GType[] ch, boolean isXlike) {
    mSeq = seq;
    mPosn = posn;
    mFa = fa;
    mMo = mo;
    mCh = ch;
    mChP = childPatterns(mFa, mMo, mCh);
    mIsXlike = isXlike;
  }

  String seq() {
    return mSeq;
  }

  int posn() {
    return mPosn;
  }

  int length() {
    return mCh.length;
  }

  Pattern pattern(final int i) {
    return mChP.index(i);
  }

  PatternArray pattern() {
    return mChP;
  }

  boolean isXLike() {
    return mIsXlike;
  }

  GType father() {
    return mFa;
  }

  GType mother() {
    return mMo;
  }

  GType child(int index) {
    return mCh[index];
  }

  /**
   * @return true iff all the children are Mendelian consistent.
   */
  boolean isMendelian() {
    for (GType aMCh : mCh) {
      final GTypeMendelian mendelianChecker = GTypeMendelianFactory.getGTypeMendelian(mFa.ploidy(), mMo.ploidy(), aMCh.ploidy());
      if (mendelianChecker == null) {
        throw new IllegalArgumentException("The ploidy combination father: " + mFa.ploidy() + " mother: " + mMo.ploidy() + " child: " + aMCh.ploidy() + " is not supported.");
      }
      if (!mendelianChecker.isMendelian(mFa, mMo, aMCh)) {
        return false;
      }
    }
    return true;
  }

  /**
   * @return true iff each parent has a singular allele (and thus can supply no phasing information).
   */
  boolean parentsSingleAllele() {
    return mFa.isSingleAllele() && mMo.isSingleAllele();
  }

  /**
   * @return true iff all the samples have ploidy none
   */
  boolean allPloidyNone() {
    return mFa.ploidy() == Ploidy.NONE && mMo.ploidy() == Ploidy.NONE;
  }

  /**
   * @return true iff both the parents and children are all heterozygous (with the same genotype). This is unlikely and can supply no phasing information.
   */
  boolean isAllHeterozygous() {
    if (mFa.ploidy() != Ploidy.DIPLOID || mMo.ploidy() != Ploidy.DIPLOID) {
      return false;
    }
    final int a = mFa.a();
    final int b = mFa.b();
    if (a == b) {
      return false;
    }
    final boolean eq0 = mMo.a() == a && mMo.b() == b;
    if (!eq0) {
      return false;
    }
    for (final GType gt : mCh) {
      final boolean eq = gt.a() == a && gt.b() == b;
      if (!eq) {
        return false;
      }
    }
    return true;
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mSeq != null && mSeq.trim().length() > 0);
    Exam.assertTrue(mPosn >= 1);
    Exam.assertNotNull(mFa);
    Exam.assertNotNull(mMo);
    Exam.assertNotNull(mCh);
    Exam.assertTrue(mCh.length > 1);
    Exam.assertEquals(mCh.length, mChP.length());
    return true;
  }

  @Override
  public void toString(StringBuilder sb) {
    sb.append(mSeq);
    sb.append(" ").append(mPosn);
    sb.append(" ").append(mFa.toString());
    sb.append(" ").append(mMo.toString());
    for (GType aMCh : mCh) {
      sb.append(" ").append(aMCh.toString());
    }
  }
}
