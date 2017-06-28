/*
 * Copyright (c) 2017. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */
package com.rtg.variant.bayes.complex;

import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import com.rtg.util.StringUtils;
import com.rtg.variant.Variant;
import com.rtg.variant.VariantLocus;
import com.rtg.variant.VariantSample;
import com.rtg.variant.bayes.AlleleStatistics;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.multisample.VariantAlleleTrigger;
import com.rtg.variant.bayes.snp.DescriptionNone;
import com.rtg.variant.util.VariantUtils;

/**
 * Some common functionality of decomposers.
 */
public abstract class AbstractDecomposer implements Decomposer {

  /**
   * Extract all the non-reference alleles from the variant record.
   * Note this is restricted to the alleles actually used in the call and need not
   * include all the alleles in the underlying description.
   * @param variant variant object
   * @return set of alleles
   */
  protected static Set<String> extractAlts(final Variant variant) {
    final Set<String> alleles = new LinkedHashSet<>();
    for (int k = 0; k < variant.getNumberOfSamples(); ++k) {
      final VariantSample vs = variant.getSample(k);
      if (vs != null) {
        if (!vs.isIdentity()) {
          Collections.addAll(alleles, StringUtils.split(vs.getName(), VariantUtils.COLON));
        }
        if (vs.getVariantAllele() != null) {
          alleles.add(vs.getVariantAllele());
        }
      }
    }
    return alleles;
  }

  /**
   * Extract all the alleles from the variant record including the reference allele.
   * The reference allele will be the first allele.
   * @param variant variant object
   * @return set of alleles
   */
  protected static String[] extractAlleles(final Variant variant) {
    final Set<String> alleles = extractAlts(variant);
    final String ref = variant.getLocus().getRefNts();
    final String[] res = new String[alleles.size() + 1];
    res[0] = ref;
    int k = 0;
    for (final String a : alleles) {
      res[++k] = a;
    }
    return res;
  }

  protected static char getAnchorBase(final VariantLocus locus, final int offset) {
    return offset == 0 ? locus.getPreviousRefNt() : locus.getRefNts().charAt(offset - 1);
  }

  static String findVariantAllele(final VariantSample sample, final String refNts, final VariantAlleleTrigger variantAlleleTrigger) {
    final AlleleStatistics<?> counts = sample.getStats().counts();
    final Description description = counts.getDescription();
    if (description instanceof DescriptionNone) {
      return null;
    }
    final int va = variantAlleleTrigger.getVariantAllele(counts, description, refNts);
    return va == -1 ? null : description.name(va);
  }

  /**
   * Changes the de novo flag from true to false on any sample incorrectly marked as de novo.
   * (Note does not do the opposite.)
   * Also deals with the somatic cause. Some variants after splitting will no longer have a somatic cause.
   * @param variant the variant to check
   * @param checker the de novo checker for the current inheritance scenario
   * @return the corrected variant
   */
  public static Variant denovoCorrect(DenovoChecker checker, Variant variant) {
    final Set<Integer> nonDenovoSamples = new HashSet<>();
    for (int s = 0; s < variant.getNumberOfSamples(); ++s) {
      final VariantSample sample = variant.getSample(s);
      if (sample != null) {
        final VariantSample.DeNovoStatus denovoCall = sample.isDeNovo();
        if (denovoCall ==  VariantSample.DeNovoStatus.IS_DE_NOVO) {
          if (!checker.isDenovo(variant, s)) {
            nonDenovoSamples.add(s);
          }
        }
      }
    }
    final Variant ret;
    if (nonDenovoSamples.size() == 0) {
      ret = variant;
    } else {
      final VariantLocus newLocus = variant.getLocus();
      final VariantSample[] newSamples = new VariantSample[variant.getNumberOfSamples()];
      for (int k = 0; k < newSamples.length; ++k) {
        final VariantSample sample = variant.getSample(k);
        if (sample != null) {
          final VariantSample.DeNovoStatus newStatus;
          if (sample.isDeNovo() == VariantSample.DeNovoStatus.UNSPECIFIED) {
            newStatus = VariantSample.DeNovoStatus.UNSPECIFIED;
          } else {
            newStatus = sample.isDeNovo() == VariantSample.DeNovoStatus.IS_DE_NOVO && !nonDenovoSamples.contains(k)
              ? VariantSample.DeNovoStatus.IS_DE_NOVO : VariantSample.DeNovoStatus.NOT_DE_NOVO;
          }
          newSamples[k] = new VariantSample(sample.getPloidy(), sample.getName(), sample.isIdentity(), sample.getMeasure(), newStatus, sample.getDeNovoPosterior());
          VariantSample.copy(sample, newSamples[k]);
        }
      }
      final Variant newVariant = new Variant(newLocus, newSamples);
      Variant.copy(variant, newVariant);
      ret = newVariant;
    }
    return ret;
  }

  protected final DenovoChecker mDenovoChecker;
  protected final VariantAlleleTrigger mVariantAlleleTrigger;

  protected AbstractDecomposer(final DenovoChecker denovoChecker, final VariantAlleleTrigger variantAlleleTrigger) {
    mDenovoChecker = denovoChecker;
    mVariantAlleleTrigger = variantAlleleTrigger;
  }

  @Override
  public abstract List<Variant> decompose(Variant original);
}
