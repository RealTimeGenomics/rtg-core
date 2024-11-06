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

package com.rtg.variant.bayes.multisample.family;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.rtg.relation.Family;
import com.rtg.util.Utils;
import com.rtg.variant.bayes.Code;
import com.rtg.variant.bayes.CodeDiploid;
import com.rtg.vcf.VcfAnnotator;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.VcfUtils;
import com.rtg.vcf.header.InfoField;
import com.rtg.vcf.header.MetaType;
import com.rtg.vcf.header.VcfHeader;
import com.rtg.vcf.header.VcfNumber;

/**
 * Class for annotating a family caller VCF output with a log segregation score.
 */
public class SegregationVcfAnnotator implements VcfAnnotator {

  /** VCF info field name for segregation probability */
  public static final String NAME = "SGP";

  private final Family mFamily;
  private final int mDecimalPlaces;

  private final Map<Integer, Code> mCodes = new HashMap<>();

  private int mFatherIndex;
  private int mMotherIndex;

  private List<Integer> mChildrenIndexes;


  /**
   * Constructor
   * @param family the family to produce a segregation number for
   */
  public SegregationVcfAnnotator(Family family) {
    this(family, 3);
  }

  /**
   * Constructor
   * @param family the family to produce a segregation number for
   * @param decimalPlaces the number of decimal places to output
   */
  public SegregationVcfAnnotator(Family family, int decimalPlaces) {
    mFamily = family;
    mDecimalPlaces = decimalPlaces;
  }

  static boolean checkHeader(VcfHeader header, Family family) {
    final int fatherIndex = header.getSampleIndex(family.getFather());
    final int motherIndex = header.getSampleIndex(family.getMother());
    boolean ret = fatherIndex != -1 && motherIndex != -1 && fatherIndex != motherIndex;
    final Set<Integer> indexes = new HashSet<>(family.size());
    indexes.add(fatherIndex);
    indexes.add(motherIndex);
    for (final String child : family.getChildren()) {
      final int childIndex = header.getSampleIndex(child);
      ret &= childIndex != -1 && !indexes.contains(childIndex);
      indexes.add(childIndex);
    }
    return ret;
  }

  @Override
  public void updateHeader(VcfHeader header) {
    assert checkHeader(header, mFamily);
    if (mChildrenIndexes == null) {
      final String[] children = mFamily.getChildren();
      mChildrenIndexes = new ArrayList<>(children.length);
      mFatherIndex = header.getSampleIndex(mFamily.getFather());
      mMotherIndex = header.getSampleIndex(mFamily.getMother());
      for (final String child : children) {
        final int childIndex = header.getSampleIndex(child);
        if (childIndex == -1) {
          throw new NullPointerException(child);
        }
        mChildrenIndexes.add(childIndex);
      }
    }
    header.ensureContains(new InfoField(NAME, MetaType.FLOAT, VcfNumber.ONE, "Segregation probability"));
  }

  @Override
  public void annotate(VcfRecord rec) {
    assert mChildrenIndexes != null && mFatherIndex >= 0 && mMotherIndex >= 0 && mFatherIndex != mMotherIndex;
    final int codeSize = rec.getAltCalls().size() + 1;
    final Code code = mCodes.computeIfAbsent(codeSize, CodeDiploid::new);
    final List<String> calls = rec.getFormat(VcfUtils.FORMAT_GENOTYPE);
    final int father = getCode(code, calls.get(mFatherIndex));
    final int mother = getCode(code, calls.get(mMotherIndex));
    if (father == -1 && mother == -1) {
      return; //Both of the parents are missing
    }

    final boolean fatherHaploid = VcfUtils.isHaploid(rec, mFatherIndex);
    final boolean motherHaploid = VcfUtils.isHaploid(rec, mMotherIndex);

    final ISegregationScore score = getSegregationScore(code, father, mother, fatherHaploid, motherHaploid);
    for (final int childIndex : mChildrenIndexes) {
      final int child = getCode(code, calls.get(childIndex));
      if (child >= 0) {
        score.increment(child, VcfUtils.isDiploid(rec, childIndex));
      }
    }
    rec.setInfo(NAME, Utils.realFormat(score.lnProbability(), mDecimalPlaces));
  }

  static ISegregationScore getSegregationScore(Code code, final int father, final int mother, final boolean fatherHaploid, final boolean motherHaploid) {
    final ISegregationScore score;
    if (father == -1 || mother == -1) {
      score = new SegregationTrivial();
    } else if (fatherHaploid) {
      if (motherHaploid) {
        score = new SegregationTrivial();
      } else {
        score = AbstractSegregationHaploid.getHaploidInstance(code, father, mother);
      }
    } else if (motherHaploid) {
      score = AbstractSegregationHaploid.getHaploidInstance(code, mother, father);
    } else {
      score = new SegregationScore(code, father, mother);
    }
    return score;
  }

  private int getCode(Code code, String call) {
    final int[] alleles = VcfUtils.splitGt(call);
    if (alleles[0] == -1) {
      return -1;
    }
    if (alleles.length == 2) {
      return code.code(alleles[0], alleles[1]);
    }
    return code.code(alleles[0]);
  }
}
