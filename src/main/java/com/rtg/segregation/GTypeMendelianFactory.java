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

import java.util.HashMap;
import java.util.Map;

import com.rtg.reference.Ploidy;
import com.rtg.util.Utils;

/**
 */
public final class GTypeMendelianFactory {

  private GTypeMendelianFactory() { }

  static final class TriplePloidy {
    final Ploidy mFa;
    final Ploidy mMo;
    final Ploidy mCh;
    TriplePloidy(Ploidy fa, Ploidy mo, Ploidy ch) {
      mFa = fa;
      mMo = mo;
      mCh = ch;
    }
    @Override
    public int hashCode() {
      return Utils.pairHashContinuous(mFa.hashCode(), mMo.hashCode(), mCh.hashCode());
    }
    @Override
    public boolean equals(Object obj) {
      if (obj == null || !(obj instanceof TriplePloidy)) {
        return false;
      }
      final TriplePloidy other = (TriplePloidy) obj;
      return other.mFa == this.mFa && other.mMo == this.mMo && other.mCh == this.mCh;
    }
  }

  private static final class PloidyMap {

    final Map<TriplePloidy, GTypeMendelian> mImplementations = new HashMap<>();

    void put(Ploidy fa, Ploidy mo, Ploidy ch, GTypeMendelian imp) {
      mImplementations.put(new TriplePloidy(fa, mo, ch), imp);
    }

    GTypeMendelian get(Ploidy fa, Ploidy mo, Ploidy ch) {
      return mImplementations.get(new TriplePloidy(fa, mo, ch));
    }

  }

  private static final PloidyMap IMPLEMENTATIONS = new PloidyMap();

  /**
   * @param fa ploidy of the father
   * @param mo ploidy of the mother
   * @param ch ploidy of the child
   * @return the GTypeMendelian implementation for the scenario
   */
  public static GTypeMendelian getGTypeMendelian(Ploidy fa, Ploidy mo, Ploidy ch) {
    return IMPLEMENTATIONS.get(fa, mo, ch);
  }

  static final class Swapped implements GTypeMendelian {
    private final GTypeMendelian mInternal;
    Swapped(GTypeMendelian internal) {
      mInternal = internal;
    }
    @Override
    public boolean isMendelian(GType fa, GType mo, GType ch) {
      return mInternal.isMendelian(mo, fa, ch);
    }
  }

  private static final GTypeMendelian DDD = new GTypeMendelian() {
    @Override
    public boolean isMendelian(GType fa, GType mo, GType ch) {
      assert fa.ploidy() == Ploidy.DIPLOID && mo.ploidy() == Ploidy.DIPLOID && ch.ploidy() == Ploidy.DIPLOID;
      if (fa.a() == ch.a()) {
        if (ch.b() == mo.a() || ch.b() == mo.b()) {
          return true;
        }
      }
      if (fa.a() == ch.b()) {
        if (ch.a() == mo.a() || ch.a() == mo.b()) {
          return true;
        }
      }
      if (fa.b() == ch.a()) {
        if (ch.b() == mo.a() || ch.b() == mo.b()) {
          return true;
        }
      }
      if (fa.b() == ch.b()) {
        if (ch.a() == mo.a() || ch.a() == mo.b()) {
          return true;
        }
      }
      return false;
    }
  };

  private static final GTypeMendelian HDD = new GTypeMendelian() {
    @Override
    public boolean isMendelian(GType fa, GType mo, GType ch) {
      assert fa.ploidy() == Ploidy.HAPLOID && mo.ploidy() == Ploidy.DIPLOID && ch.ploidy() == Ploidy.DIPLOID;
      if (fa.a() == ch.a()) {
        if (ch.b() == mo.a() || ch.b() == mo.b()) {
          return true;
        }
      }
      if (fa.a() == ch.b()) {
        if (ch.a() == mo.a() || ch.a() == mo.b()) {
          return true;
        }
      }
      return false;
    }
  };

  private static final GTypeMendelian HDH = new GTypeMendelian() {
    @Override
    public boolean isMendelian(GType fa, GType mo, GType ch) {
      assert fa.ploidy() == Ploidy.HAPLOID && mo.ploidy() == Ploidy.DIPLOID && ch.ploidy() == Ploidy.HAPLOID;
      return ch.a() == mo.a() || ch.a() == mo.b();
    }
  };

  private static final GTypeMendelian NHN = new GTypeMendelian() {
    @Override
    public boolean isMendelian(GType fa, GType mo, GType ch) {
      assert fa.ploidy() == Ploidy.NONE && mo.ploidy() == Ploidy.HAPLOID && ch.ploidy() == Ploidy.NONE;
      return true;
    }
  };

  private static final GTypeMendelian HNH = new GTypeMendelian() {
    @Override
    public boolean isMendelian(GType fa, GType mo, GType ch) {
      assert fa.ploidy() == Ploidy.HAPLOID && mo.ploidy() == Ploidy.NONE && ch.ploidy() == Ploidy.HAPLOID;
      return fa.a() == ch.a();
    }
  };

  private static final GTypeMendelian PPP = new GTypeMendelian() {
    @Override
    public boolean isMendelian(GType fa, GType mo, GType ch) {
      assert fa.ploidy() == Ploidy.POLYPLOID && mo.ploidy() == Ploidy.POLYPLOID && ch.ploidy() == Ploidy.POLYPLOID;
      return mo.a() == ch.a();
    }
  };

  static {
    final Ploidy dip = Ploidy.DIPLOID;
    final Ploidy non = Ploidy.NONE;
    final Ploidy hap = Ploidy.HAPLOID;
    final Ploidy pol = Ploidy.POLYPLOID;
    IMPLEMENTATIONS.put(dip, dip, dip, DDD);
    IMPLEMENTATIONS.put(hap, dip, dip, HDD);
    IMPLEMENTATIONS.put(dip, hap, dip, new Swapped(HDD));
    IMPLEMENTATIONS.put(hap, dip, hap, HDH);
    IMPLEMENTATIONS.put(dip, hap, hap, new Swapped(HDH));
    IMPLEMENTATIONS.put(non, hap, non, NHN);
    IMPLEMENTATIONS.put(hap, non, non, new Swapped(NHN));
    IMPLEMENTATIONS.put(hap, non, hap, HNH);
    IMPLEMENTATIONS.put(non, hap, hap, new Swapped(HNH));
    IMPLEMENTATIONS.put(pol, pol, pol, PPP);
  }
}
