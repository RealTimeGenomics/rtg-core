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
package com.rtg.variant.bayes.complex;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.rtg.util.StringUtils;
import com.rtg.variant.Variant;
import com.rtg.variant.VariantLocus;
import com.rtg.variant.VariantOutputLevel;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantSample;
import com.rtg.variant.bayes.AlleleStatistics;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Statistics;
import com.rtg.variant.bayes.multisample.VariantAlleleTrigger;
import com.rtg.variant.bayes.snp.DescriptionCommon;
import com.rtg.variant.bayes.snp.DescriptionNone;
import com.rtg.variant.util.VariantUtils;
import com.rtg.variant.util.arithmetic.LogApproximatePossibility;

/**
 * Set of utilities to deal with trimming of hypotheses from the complex caller.
 */
public final class Trimming {

  private Trimming() { }

  // Can't use Description, as that includes hypotheses that weren't necessarily called.
  private static HashSet<String> extractAlleles(final Variant variant) {
    final HashSet<String> alleles = new HashSet<>();
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

  private static VariantLocus createLocus(final Variant original, final int leftClip, final int rightClip) {
    final VariantLocus locus = original.getLocus();
    final char newPrevNt = leftClip == 0 ? locus.getPreviousRefNt() : locus.getRefNts().charAt(leftClip - 1);
    final String newReference = StringUtils.clip(locus.getRefNts(), leftClip, rightClip);
    final int newStart = locus.getStart() + leftClip;
    final int newEnd = locus.getEnd() - rightClip;
    return new VariantLocus(locus.getSequenceName(), newStart, newEnd, newReference, newPrevNt);
  }

  private static String createVariantName(final int leftClip, final int rightClip, final String oldName) {
    if (oldName == null) {
      return oldName;
    }
    final int colon = oldName.indexOf(VariantUtils.COLON);
    if (colon >= 0) {
      assert oldName.indexOf(VariantUtils.COLON, colon + 1) == -1;
      final String allele1 = StringUtils.clip(oldName.substring(0, colon), leftClip, rightClip);
      final String allele2 = StringUtils.clip(oldName.substring(colon + 1), leftClip, rightClip);
      return allele1 + VariantUtils.COLON + allele2;
    } else {
      assert oldName.indexOf(VariantUtils.COLON) == -1;
      return StringUtils.clip(oldName, leftClip, rightClip);
    }
  }

  private static String findVariantAllele(VariantSample sample, String refNts, VariantAlleleTrigger variantAlleleTrigger) {
    final AlleleStatistics<?> counts = sample.getStats().counts();
    final Description description = counts.getDescription();
    if (description instanceof DescriptionNone) {
      return null;
    }
    final int va = variantAlleleTrigger.getVariantAllele(counts, description, refNts);
    return va == -1 ? null : description.name(va);
  }

  private static VariantSample[] createVariants(final Variant original, final int leftClip, final int rightClip, Description newDescription, int[] alleleMapping, String refNts, VariantAlleleTrigger variantAlleleTrigger) {
    final VariantSample[] newSamples = new VariantSample[original.getNumberOfSamples()];
    for (int k = 0; k < newSamples.length; ++k) {
      final VariantSample sample = original.getSample(k);
      if (sample != null) {
        final String newName = createVariantName(leftClip, rightClip, sample.getName());
        newSamples[k] = new VariantSample(sample.getPloidy(), newName, sample.isIdentity(), sample.getMeasure(), sample.isDeNovo(), sample.getDeNovoPosterior());
        VariantSample.copy(sample, newSamples[k]);
        // Update counts to correspond to the new description
        final Description d = newSamples[k].getStats().counts().getDescription() instanceof DescriptionNone ? DescriptionNone.SINGLETON : newDescription;
        newSamples[k].setStats((Statistics<?>) newSamples[k].getStats().copy());
        newSamples[k].getStats().remapAlleleStatistics(d, alleleMapping);
        newSamples[k].setVariantAllele(findVariantAllele(newSamples[k], refNts, variantAlleleTrigger));
        final Map<Set<String>, Double> newMap = newGenotypeLikelihoods(leftClip, rightClip, sample);
        newSamples[k].setGenotypeLikelihoods(newMap);
      } else {
        newSamples[k] = null;
      }
    }
    return newSamples;
  }

  // XXX These should be recomputed directly from the measure, in order to include alleles not called in the originalLikelihoods
  private static Map<Set<String>, Double> newGenotypeLikelihoods(int leftClip, int rightClip, VariantSample sample) {
    final Map<Set<String>, Double> newMap = new HashMap<>();
    final Map<Set<String>, Double> originalLikelihoods = sample.getGenotypeLikelihoods();
    if (originalLikelihoods != null) {
      for (Map.Entry<Set<String>, Double> entry : originalLikelihoods.entrySet()) {
        final Set<String> newSet = new HashSet<>();
        for (String s : entry.getKey()) {
          newSet.add(StringUtils.clip(s, leftClip, rightClip));
        }
        final Double v = newMap.get(newSet);
        final double existing = v == null ? LogApproximatePossibility.SINGLETON.zero() : v;
        newMap.put(newSet, LogApproximatePossibility.SINGLETON.add(existing, entry.getValue()));
      }
      return newMap;
    } else {
      return null;
    }
  }


  /**
   * Trim a variant by removing common prefix and suffix from each call and reference
   * and adjusting all the position information accordingly.
   * @param original the original call
   * @param variantAlleleTrigger variant allele detector
   * @return trimmed variant (possibly original if no change made)
   */
  static Variant trim(final Variant original, VariantAlleleTrigger variantAlleleTrigger) {

    final String ref = original.getLocus().getRefNts();
    if (ref.length() == 0) {
      // Cannot possibly trim if we are inserting, all bases are to be inserted
      return original;
    }

    // Compute set of called alleles and exit if there are none
    final HashSet<String> catSet = extractAlleles(original);
    if (catSet.size() == 0) {
      return original;
    }

    // Include actual reference sequence (this will not be "")
    assert ref.length() > 0;
    catSet.add(ref);

    // Compute maximal clip positions based on set of alleles
    final String[] cats = catSet.toArray(new String[catSet.size()]);

    final int rightClip = StringUtils.longestSuffix(cats);
    final int leftClip = StringUtils.longestPrefix(rightClip, cats);
    //final int leftClip = StringUtils.longestPrefix(cats);
    //final int rightClip = StringUtils.longestSuffix(cats, leftClip);

    // Quick exit if no trimming is possible
    if (leftClip == 0 && rightClip == 0) {
      return original;
    }

    // Create new locus, will be "" for case of pure insertion
    final VariantLocus newLocus = createLocus(original, leftClip, rightClip);

    final Variant result = getVariant(original, leftClip, rightClip, variantAlleleTrigger, newLocus);
    result.setTrimmed();
    //System.err.println(original);
    //System.err.println(result);
    return result;
  }

  private static Variant createSplitVariant(final Variant original, final int start, final int end, final int id, VariantAlleleTrigger variantAlleleTrigger) {
    final VariantLocus newLocus = createLocus(original, start, end);
    final Variant result = getVariant(original, start, end, variantAlleleTrigger, newLocus);
    result.setSplitId(id);
    return result;
  }

  private static Variant getVariant(Variant original, int start, int end, VariantAlleleTrigger variantAlleleTrigger, VariantLocus newLocus) {
    final VariantSample[] newSamples1;
    if (original.getNumberOfSamples() > 0) {
      //trim description
      Description oldDescription = DescriptionNone.SINGLETON;
      for (int i = 0; i < original.getNumberOfSamples(); ++i) {
        if (original.getSample(i) != null && !(original.getSample(i).getStats().counts().getDescription() instanceof DescriptionNone)) {
          oldDescription = original.getSample(i).getStats().counts().getDescription();
        }
      }
      // Incrementally build up mapping of old alleles to new alleles
      final LinkedHashMap<String, Integer> alleles = new LinkedHashMap<>();
      final int[] alleleMap = new int[oldDescription.size()];
      for (int i = 0; i < oldDescription.size(); ++i) {
        final String clipped = StringUtils.clip(oldDescription.name(i), start, end);
        Integer newPos = alleles.get(clipped);
        if (newPos == null) {
          newPos = alleles.size();
        }
        alleles.put(clipped, newPos);
        alleleMap[i] = newPos;
      }

      final Description newDescription = new DescriptionCommon(alleles.keySet().toArray(new String[alleles.size()]));
      newSamples1 = createVariants(original, start, end, newDescription, alleleMap, newLocus.getRefNts(), variantAlleleTrigger);
    } else {
      newSamples1 = new VariantSample[0];
    }
    final VariantSample[] newSamples = newSamples1;
    final Variant result = new Variant(newLocus, newSamples);
    Variant.copy(original, result);
    result.setPossibleCause(createVariantName(start, end, original.getPossibleCause()));
    return result;
  }

  static List<Variant> split(final Variant original, DenovoChecker denovoCorrector, VariantAlleleTrigger variantAlleleTrigger) {

    // Compute set of alleles and exit if there are none
    final HashSet<String> catSet = extractAlleles(original);
    if (catSet.size() == 0) {
      return Collections.singletonList(original);
    }
    final String ref = original.getLocus().getRefNts();
    catSet.add(ref);

    // Check all are the same length, if not we cannot split
    final String[] cats = catSet.toArray(new String[catSet.size()]);
    final int length = cats[0].length();
    for (final String c : cats) {
      if (c.length() != length) {
        return Collections.singletonList(original);
      }
    }

    // After this loop, "false" in syndrome indicates a column where all alleles agree
    // and thus represents a split point.
    final boolean[] syndrome = new boolean[length];
    final String c = cats[0];
    for (int k = 1; k < cats.length; ++k) {
      for (int j = 0; j < length; ++j) {
        syndrome[j] |= c.charAt(j) != cats[k].charAt(j);
      }
    }
    //System.err.println(java.util.Arrays.toString(syndrome));

    // Check for at least one split point
    boolean hasNoSplitPoint = true;
    for (final boolean s : syndrome) {
      hasNoSplitPoint &= s;
    }
    if (hasNoSplitPoint) {
      return Collections.singletonList(original);
    }

    // Create variants around the split points
    final List<Variant> list = new ArrayList<>();
    int startSplit = 0;
    int splitId = 0; // Unique identifier for each subcall
    while (startSplit < length) {
      if (syndrome[startSplit]) {
        int endSplit = startSplit + 1;
        while (endSplit < length && syndrome[endSplit]) {
          ++endSplit;
        }
        final Variant splitVariant = createSplitVariant(original, startSplit, length - endSplit, splitId++, variantAlleleTrigger);
        final Variant variant = denovoCorrector != null ? denovoCorrect(denovoCorrector, splitVariant) : splitVariant;
        list.add(variant);
        startSplit = endSplit + 1;
      } else {
        ++startSplit;
      }
    }
    return list;
  }

  /**
   * changes the de novo flag from true to false on any sample incorrectly marked as de novo. (note does not do the opposite)
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

  /**
   * Trims and splits complex calls.
   *
   * @param params params to check options
   * @param original call to be trimmed and split (may be returned as the answer)
   * @param denovoCorrector method for correcting de novo flags
   * @return a list of calls. There can be more than one if the call can be split into separate calls that line up.
   */
  public static List<Variant> trimSplit(VariantParams params, Variant original, DenovoChecker denovoCorrector) {
    final VariantAlleleTrigger variantAlleleTrigger = new VariantAlleleTrigger(params.minVariantAllelicDepth(), params.minVariantAllelicFraction());
    if (params.callLevel() == VariantOutputLevel.ALL) {
      return Collections.singletonList(original);
    } else if (params.ionTorrent()) {
      return Collections.singletonList(trim(original, variantAlleleTrigger));
    }
    return split(trim(original, variantAlleleTrigger), denovoCorrector, variantAlleleTrigger);
  }
}
