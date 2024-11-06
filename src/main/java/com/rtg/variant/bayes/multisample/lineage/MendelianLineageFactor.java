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
package com.rtg.variant.bayes.multisample.lineage;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.rtg.reference.Ploidy;
import com.rtg.variant.bayes.Code;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.multisample.family.UniqueId;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 * Mendelian for genotype inheritance in a cell lineage.
 */
public class MendelianLineageFactor extends AbstractFactor implements ToDefaultFactor {

  private final List<Variable> mScope = new ArrayList<>();
  private final Mendelian mMendelian;
  private final Variable mGenotype;
  private final Variable mGenotypeParent;
  private final Variable mDeNovo;

  private interface Mendelian {
    double p(int gChild, int gParent, int deNovo);
  }

  private static Mendelian haploidMendelian(final int m, final double deNovo) {
    return new Mendelian() {
      private final double[][] mP = {{1 - deNovo, 0}, {0, deNovo / (m - 1)}};

      @Override
      public double p(int gChild, int gParent, int deNovo) {
        return mP[deNovo][gChild == gParent ? 0 : 1];
      }
    };
  }

  private static Mendelian diploidMendelian(final Code code, final int m, final double deNovo) {
    return new Mendelian() {
      private final double mMu = 1 - Math.sqrt(1 - deNovo); //0.5 * deNovo;
      private final Code mCode = code;

      // See Table 2: Approximate Diploid Mutation in Patent
      @Override
      public double p(int gChild, int gParent, int deNovo) {
        final UniqueId uid = new UniqueId(mCode.rangeSize());
        final int gPa = uid.addId(mCode.a(gParent));  // 0
        final int gPb = uid.addId(mCode.bc(gParent)); // 0 or 1
        final int gCa = uid.addId(mCode.a(gChild));
        final int gCb = uid.addId(mCode.bc(gChild));
        if (gPa == gPb) {
          // parent is homozygous
          if (gCa != gPa && gCb != gPa) {
            return 0;
          }
          if (gCa == gCb /* && gPa == gCa */) {
            return deNovo == 0 ? (1 - mMu) / (1 + mMu) : 0;
          }
          // Child has one allele different to parent
          return deNovo == 0 ? 0 : 2 * mMu / ((m - 1) * (1 + mMu));
        } else {
          if (deNovo == 0) {
            if ((gCa == gPa && gCb == gPb) || (gCa == gPb && gCb == gPa)) {
              return (1 - mMu) / (1 + mMu);
            } else {
              return 0;
            }
          } else {
            // De novo is true
            // Want cases where exactly one allele in the child is different from parent
            if ((gCa == gPa && gCb != gPb) || (gCa == gPb && gCb != gPa) || (gCb == gPa && gCa != gPb) || (gCb == gPb && gCa != gPa)) {
              return mMu / ((m - 1) * (1 + mMu));
            } else {
              return 0;
            }
          }
        }
      }
    };
  }

  /**
   * @param arithmetic the arithmetic
   * @param genotype genotype variable of the child
   * @param genotypeParent genotype variable of the parent
   * @param deNovo de novo variable of the child (as a probability)
   * @param deNovoPrior the prior for a de novo mutation in this child
   * @param hypotheses genotype hypotheses for the child
   * @param hypothesesParent genotype hypotheses for the parent
   */
  public MendelianLineageFactor(PossibilityArithmetic arithmetic, Variable genotype, Variable genotypeParent, Variable deNovo, double deNovoPrior, Hypotheses<?> hypotheses, Hypotheses<?> hypothesesParent) {
    super(arithmetic);
    mGenotype = genotype;
    mGenotypeParent = genotypeParent;
    mDeNovo = deNovo;
    // This order is important for asDefault to work
    mScope.add(genotype);
    mScope.add(genotypeParent);
    mScope.add(deNovo);
    if (hypotheses.ploidy() != hypothesesParent.ploidy()) {
      throw new UnsupportedOperationException("Different ploidy in lineage not supported");
    }
    if (hypotheses.ploidy() == Ploidy.HAPLOID) {
      mMendelian = haploidMendelian(hypotheses.size(), deNovoPrior);
    } else if (hypotheses.ploidy() == Ploidy.DIPLOID) {
      assert hypotheses.code().equals(hypothesesParent.code());
      mMendelian = diploidMendelian(hypotheses.code(), hypotheses.size(), deNovoPrior);
    } else {
      throw new UnsupportedOperationException("Can't do ploidy " + hypotheses.ploidy());
    }
  }

  @Override
  public double p(Map<Variable, Integer> values) {
    if (!values.keySet().equals(scope())) {
      throw new IllegalArgumentException("Scope errror");
    }
    return arithmetic().prob2Poss(mMendelian.p(values.get(mGenotype), values.get(mGenotypeParent), values.get(mDeNovo)));
  }

  @Override
  public Set<Variable> scope() {
    return new HashSet<>(mScope);
  }

  @Override
  public DefaultFactor asDefault() {
    int factorSize = 1;
    for (final Variable s : mScope) {
      factorSize *= s.size();
    }
    final double[] poss = new double[factorSize];
    for (int deNovo = 0, k = 0; deNovo < mDeNovo.size(); ++deNovo) {
      for (int genotypeParent = 0; genotypeParent < mGenotypeParent.size(); ++genotypeParent) {
        for (int genotype = 0; genotype < mGenotype.size(); ++genotype, ++k) {
          poss[k] = arithmetic().prob2Poss(mMendelian.p(genotype, genotypeParent, deNovo));
        }
      }
    }

    return new DefaultFactor(arithmetic(), mScope, poss);
  }

}
