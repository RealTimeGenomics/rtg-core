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

package com.rtg.segregation;

import java.util.Arrays;

import com.rtg.reference.Ploidy;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 */
class Labelling extends IntegralAbstract {
  private static final int NUMBER_LABELS = 4;
  static final int MOTHER_LABEL = 2;
  static final int FATHER_LABEL = 1;

  static boolean equivalent(final GType x, final GType y) {
    if (x.a() == y.a() && x.b() == y.b()) {
      return true;
    }
    if (x.a() == y.b() && x.b() == y.a()) {
      return true;
    }
    return false;
  }

  static boolean aSame(final GType x, final GType y) {
    if (x.a() == y.a()) {
      return true;
    }
    if (x.a() == y.b()) {
      return true;
    }
    return false;
  }

  static Integer label(final PatternArray pattern, final int index) {
    final String fa = pattern.faString().substring(index, index + 1);
    final int l0;
    switch(fa) {
      case "0":
        l0 = 0;
        break;
      case "1":
        l0 = FATHER_LABEL;
        break;
      default:
        return null;
    }
    final String mo = pattern.moString().substring(index, index + 1);
    final int l1;
    switch(mo) {
      case "0":
        l1 = 0;
        break;
      case "1":
        l1 = MOTHER_LABEL;
        break;
      default:
        return null;
    }
    return l0 + l1;
  }

  /**
   * Phase without using labels if possible.
   * @param father genotype.
   * @param mother genotype.
   * @param child genotype.
   * @return phased genotype for child or null if cant do it.
   */
  static GType simplePhasing(GType father, GType mother, GType child) {
    assert !child.isSingleAllele();
    final int a = child.a();
    final int b = child.b();
    final GType childAlt = new GType(b, a, child.ploidy());
    final GType aGt = simpleCase(a, father, mother, child, childAlt);
    if (aGt != null) {
      return aGt;
    }
    final GType bGt = simpleCase(b, father, mother, childAlt, child);
    if (bGt != null) {
      return bGt;
    }
    return null;
  }

  static GType simpleCase(final int x, final GType father, final GType mother, final GType primary, final GType secondary) {
    final boolean inFa = isIn(x, father);
    final boolean inMo = isIn(x, mother);
    if (inFa && !inMo) {
      return primary;
    }
    if (inMo && !inFa) {
      return secondary;
    }
    return null;
  }

  static boolean isIn(final int x, final GType gt) {
    return x == gt.a() || x == gt.b();
  }

  private final FamilyGt mFamily;
  private final PatternArray mPattern;
  private final GType[] mLabels;
  private final Integer[] mChildIndex;
  private final GType[] mPhasedLabels;

  Labelling(FamilyGt family, PatternArray pattern) {
    mFamily = family;
    mPattern = pattern;
    mChildIndex = new Integer[NUMBER_LABELS];
    mPhasedLabels = new GType[NUMBER_LABELS];
    mLabels = new GType[NUMBER_LABELS];
    for (int i = 0; i < family.length(); ++i) {
      final Integer label = label(pattern, i);
      if (label == null) {
        continue;
      }
      final GType child = family.child(i);
      if (mLabels[label] != null) {
        assert equivalent(mLabels[label], child);
      } else {
        mLabels[label] = child;
        mChildIndex[label] = i;
      }
    }
    //this needs to be done after the mLabels are constructed
    for (int i = 0; i < NUMBER_LABELS; ++i) {
      final GType label = mLabels[i];
      if (label != null) {
        if (label.isSingleAllele()) {
          mPhasedLabels[i] = label;
        } else {
          mPhasedLabels[i] = phaseChild(mChildIndex[i]);
        }
      }
    }
  }

  /**
   * @return the correctly phased genotype of father (0 label first, 1 label second). null if phasing not possible.
   */
  public GType phaseFather() {
    final GType father = mFamily.father();
    if (father.isSingleAllele()) {
      return null;
    }
    final GType fatherAlt = new GType(father.b(), father.a(), father.ploidy());
    for (int i = 0; i < NUMBER_LABELS; ++i) {
      final GType label = mPhasedLabels[i];
      if (label == null || label.ploidy() != Ploidy.DIPLOID) {
        continue;
      }
      final int fal = i & FATHER_LABEL;
      final int a = label.a();
      if (father.a() == a) {
        return fal == 0 ? father : fatherAlt;
      } else if (father.b() == a) {
        return fal == 0 ? fatherAlt : father;
      } else {
        throw new RuntimeException("i:" + i + " father:" + father + " phasedLabels:" + Arrays.toString(mPhasedLabels));
      }
    }
    return null;
  }

  /**
   * @return the correctly phased genotype of mother (0 label first, 1 label second). null if phasing not possible.
   */
  public GType phaseMother() {
    final GType mother = mFamily.mother();
    if (mother.isSingleAllele()) {
      return null;
    }
    final GType motherAlt = new GType(mother.b(), mother.a(), mother.ploidy());
    for (int i = 0; i < NUMBER_LABELS; ++i) {
      final GType label = mPhasedLabels[i];
      if (label == null || label.ploidy() != Ploidy.DIPLOID) {
        continue;
      }
      final int mol = i & MOTHER_LABEL;
      final int b = label.b();
      if (mother.a() == b) {
        return mol == 0 ? mother : motherAlt;
      } else if (mother.b() == b) {
        return mol == 0 ? motherAlt : mother;
      } else {
        throw new RuntimeException("i:" + i + " mother:" + mother + " phasedLabels:" + Arrays.toString(mPhasedLabels));
      }
    }
    return null;
  }

  /**
   * @param childIndex index of child (0 based).
   * @return the correctly phased genotype of child (father first, mother second). null if phasing not possible.
   */
  public final GType phaseChild(int childIndex) {
    final GType child = mFamily.child(childIndex);
    final GType simple = simplePhasing(mFamily.father(), mFamily.mother(), child);
    final GType res;
    if (simple != null) {
      res = simple;
    } else {
      res = labelPhasing(childIndex);
    }
    return res;
  }

  GType labelPhasing(final int childIndex) {
    final GType child  = mFamily.child(childIndex);
    assert !child.isSingleAllele();
    final Integer label = label(mPattern, childIndex);
    if (label == null) {
      return null;
    }
    final GType res;
    final GType childPh = new GType(child.b(), child.a(), child.ploidy());
    final GType faGt = mLabels[label ^ FATHER_LABEL];
    if (faGt != null && !equivalent(child, faGt)) {
      res = aSame(child, faGt) ? childPh : child;
    } else {
      final GType moGt = mLabels[label ^ MOTHER_LABEL];
      if (moGt != null && !equivalent(child, moGt)) {
        res = aSame(child, moGt) ? child : childPh;
      } else {
        res = null;
      }
    }
    return res;
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    for (int i = 0; i < NUMBER_LABELS; ++i) {
      final GType label = mLabels[i];
      if (label == null) {
        Exam.assertTrue(mPhasedLabels[i] == null);
        Exam.assertTrue(mChildIndex[i] == null);
      } else {
        Exam.assertTrue(equivalent(label, mPhasedLabels[i]));
        Exam.assertTrue(label == mFamily.child(mChildIndex[i]));
        Exam.assertTrue(0 <= mChildIndex[i] && mChildIndex[i] < mFamily.length());
      }
    }
    return true;
  }

  @Override
  public boolean integrity() {
    Exam.assertNotNull(mFamily);
    Exam.assertNotNull(mPattern);
    Exam.assertEquals(mFamily.length(), mPattern.length());
    Exam.assertEquals(NUMBER_LABELS, mLabels.length);
    Exam.assertEquals(NUMBER_LABELS, mPhasedLabels.length);
    Exam.assertEquals(NUMBER_LABELS, mChildIndex.length);
    return true;
  }

}
