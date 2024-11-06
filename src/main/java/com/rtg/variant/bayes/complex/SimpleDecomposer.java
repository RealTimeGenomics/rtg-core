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
import com.rtg.variant.VariantSample;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.multisample.VariantAlleleTrigger;
import com.rtg.variant.bayes.snp.DescriptionCommon;
import com.rtg.variant.bayes.snp.DescriptionNone;
import com.rtg.variant.util.VariantUtils;
import com.rtg.variant.util.arithmetic.LogApproximatePossibility;

/**
 * Set of utilities to deal with trimming of hypotheses from the complex caller.
 */
public class SimpleDecomposer extends AbstractDecomposer {

  private static VariantLocus createLocus(final Variant original, final int leftClip, final int rightClip) {
    final VariantLocus locus = original.getLocus();
    final String newReference = StringUtils.clip(locus.getRefNts(), leftClip, rightClip);
    final int newStart = locus.getStart() + leftClip;
    final int newEnd = locus.getEnd() - rightClip;
    return new VariantLocus(locus.getSequenceName(), newStart, newEnd, newReference, getAnchorBase(locus, leftClip));
  }

  private static String createVariantGenotype(final int leftClip, final int rightClip, final String oldGenotype) {
    if (oldGenotype == null) {
      return oldGenotype;
    }
    final int colon = oldGenotype.indexOf(VariantUtils.COLON);
    if (colon >= 0) {
      assert oldGenotype.indexOf(VariantUtils.COLON, colon + 1) == -1;
      final String allele1 = StringUtils.clip(oldGenotype.substring(0, colon), leftClip, rightClip);
      final String allele2 = StringUtils.clip(oldGenotype.substring(colon + 1), leftClip, rightClip);
      return allele1 + VariantUtils.COLON + allele2;
    } else {
      return StringUtils.clip(oldGenotype, leftClip, rightClip);
    }
  }

  // TODO These should be recomputed directly from the measure, in order to include alleles not called in the originalLikelihoods
  private static Map<Set<String>, Double> newGenotypeLikelihoods(final VariantSample sample, final int leftClip, final int rightClip) {
    final Map<Set<String>, Double> originalLikelihoods = sample.getGenotypeLikelihoods();
    if (originalLikelihoods != null) {
      final Map<Set<String>, Double> newMap = new HashMap<>(originalLikelihoods.size());
      for (final Map.Entry<Set<String>, Double> entry : originalLikelihoods.entrySet()) {
        final Set<String> key = entry.getKey();
        final Set<String> newSet = new HashSet<>(key.size());
        for (final String s : key) {
          newSet.add(StringUtils.clip(s, leftClip, rightClip)); // TODO assumptions about all alleles same length ??
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

  private static VariantSample[] createVariants(final Variant original, final int leftClip, final int rightClip, Description newDescription, int[] alleleMapping, String refNts, VariantAlleleTrigger variantAlleleTrigger) {
    final VariantSample[] newSamples = new VariantSample[original.getNumberOfSamples()];
    for (int k = 0; k < newSamples.length; ++k) {
      final VariantSample sample = original.getSample(k);
      if (sample != null) {
        final String newName = createVariantGenotype(leftClip, rightClip, sample.getName());
        newSamples[k] = new VariantSample(sample.getPloidy(), newName, sample.isIdentity(), sample.getMeasure(), sample.isDeNovo(), sample.getDeNovoPosterior());
        VariantSample.copy(sample, newSamples[k]);
        // Update counts to correspond to the new description
        final Description d = newSamples[k].getStats().counts().getDescription() instanceof DescriptionNone ? DescriptionNone.SINGLETON : newDescription;
        newSamples[k].getStats().remapAlleleStatistics(d, alleleMapping);
        newSamples[k].setVariantAllele(findVariantAllele(newSamples[k], refNts, variantAlleleTrigger));
        newSamples[k].setGenotypeLikelihoods(newGenotypeLikelihoods(sample, leftClip, rightClip));
      } else {
        newSamples[k] = null;
      }
    }
    return newSamples;
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
    final Set<String> catSet = extractAlts(original);
    if (catSet.isEmpty()) {
      return original;
    }

    // Include actual reference sequence (this will not be "")
    assert ref.length() > 0;
    catSet.add(ref);

    // Compute maximal clip positions based on set of alleles
    final String[] cats = catSet.toArray(new String[0]);

    final int rightClip = StringUtils.longestSuffix(cats);
    final int leftClip = StringUtils.longestPrefix(rightClip, cats);

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
    final VariantSample[] newSamples;
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
      final Description newDescription = new DescriptionCommon(alleles.keySet().toArray(new String[0]));
      //System.out.println(oldDescription + " " + newDescription + " " + Arrays.toString(alleleMap));
      newSamples = createVariants(original, start, end, newDescription, alleleMap, newLocus.getRefNts(), variantAlleleTrigger);
    } else {
      newSamples = new VariantSample[0];
    }
    final Variant result = new Variant(newLocus, newSamples);
    Variant.copy(original, result);
    result.setPossibleCause(createVariantGenotype(start, end, original.getPossibleCause()));
    return result;
  }

  /**
   * Construct a new trimmer and splitter
   * @param denovoChecker method for correcting de novo flags
   * @param variantAlleleTrigger method for handling variant alleles
   */
  public SimpleDecomposer(DenovoChecker denovoChecker, VariantAlleleTrigger variantAlleleTrigger) {
    super(denovoChecker, variantAlleleTrigger);
  }

  SimpleDecomposer() {
    this(null, new VariantAlleleTrigger(0, 0));
  }

  List<Variant> split(final Variant original, VariantAlleleTrigger variantAlleleTrigger) {

    // Compute set of alleles and exit if there are none
    final Set<String> catSet = extractAlts(original);
    if (catSet.isEmpty()) {
      return Collections.singletonList(original);
    }
    final String ref = original.getLocus().getRefNts();
    catSet.add(ref);

    // Check all are the same length, if not we cannot split
    final String[] cats = catSet.toArray(new String[0]);
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
        final Variant variant = mDenovoChecker != null ? denovoCorrect(mDenovoChecker, splitVariant) : splitVariant;
        list.add(variant);
        startSplit = endSplit + 1;
      } else {
        ++startSplit;
      }
    }
    return list;
  }

  @Override
  public List<Variant> decompose(final Variant original) {
    final List<Variant> res = split(original, mVariantAlleleTrigger);
    return original == res.get(0) ? Collections.singletonList(trim(original, mVariantAlleleTrigger)) : res;
  }
}
