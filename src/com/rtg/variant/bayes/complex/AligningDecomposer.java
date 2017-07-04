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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.rtg.alignment.SplitAlleles;
import com.rtg.util.Pair;
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
 * A decomposer that works by aligning the alternative alleles to the reference
 * during the determination of split points.
 *
 * The implementation computes a pairwise alignment of each description allele
 * to the reference.  All the alleles/cigars used in the call are then stacked
 * up and a partitioning performed using {@code SplitAlleles}.
 *
 * In order to compute better statistics, other alleles occurring in the
 * description but not the call, are partitioned separately one at a time
 * against the same reference and alternates.  This information can then be
 * used to assign the counts arising from such alleles to the split alleles of
 * the call whenever there is a match in the nucleotide sequence.
 */
public class AligningDecomposer extends AbstractDecomposer {

  /**
   * Construct a new decomposer using alignment and the specified de novo checker.
   * @param checker checker for de novo events
   * @param variantAlleleTrigger method for handling variant alleles
   */
  public AligningDecomposer(final DenovoChecker checker, final VariantAlleleTrigger variantAlleleTrigger) {
    super(checker, variantAlleleTrigger);
  }

  AligningDecomposer() {
    this(null, new VariantAlleleTrigger(0, 0));
  }

  // XXX These should be recomputed directly from the measure, in order to include alleles not called in the originalLikelihoods
  private static Map<Set<String>, Double> newGenotypeLikelihoods(final VariantSample sample, final SplitAlleles splitter, final String[] alleles, final int leftClip) {
    final Map<Set<String>, Double> newMap = new HashMap<>();
    final Map<Set<String>, Double> originalLikelihoods = sample.getGenotypeLikelihoods();
    if (originalLikelihoods != null) {
      for (final Map.Entry<Set<String>, Double> entry : originalLikelihoods.entrySet()) {
        final Set<String> newSet = new HashSet<>();
        for (final String s : entry.getKey()) {
          newSet.add(alleles[splitter.getColumnIndex(s)].substring(leftClip));
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

  private static String createVariantGenotype(final SplitAlleles splitter, final String[] part, final String oldGenotype, final int leftClip) {
    if (oldGenotype == null) {
      return oldGenotype;
    }
    final int colon = oldGenotype.indexOf(VariantUtils.COLON);
    if (colon >= 0) {
      final String allele1 = part[splitter.getColumnIndex(oldGenotype.substring(0, colon))].substring(leftClip);
      final String allele2 = part[splitter.getColumnIndex(oldGenotype.substring(colon + 1))].substring(leftClip);
      return allele1 + VariantUtils.COLON + allele2;
    } else {
      return part[splitter.getColumnIndex(oldGenotype)].substring(leftClip);
    }
  }

  private int findIfPresent(final String[] alleles, final String a) {
    for (int j = 0; j < alleles.length; ++j) {
      if (a.equals(alleles[j])) {
        return j;
      }
    }
    return -1;
  }

  private Pair<String[], int[]> getAlleleMap(final Description oldDescription, final String[] originalAlleles, final String[] extraDescriptionAlleles, final List<List<Pair<Integer, String[]>>> partitions, final int slice, final int leftClip) {
    final Pair<Integer, String[]> part = partitions.get(0).get(slice);
    final String[] alleles = part.getB();
    assert originalAlleles.length == alleles.length;
    // Incrementally build up mapping of old alleles to new alleles
    final int[] alleleMap = new int[oldDescription.size()]; // Mapping of old allele numbering onto new unique allele numbering
    Arrays.fill(alleleMap, -1);
    final LinkedHashMap<String, Integer> uniqueAlleles = new LinkedHashMap<>();
    // Deal with description alleles occurring in original call
    for (int k = 0; k < oldDescription.size(); ++k) {
      final String a = oldDescription.name(k);
      final int j = findIfPresent(originalAlleles, a);
      if (j >= 0) {
        // This description allele is reference or one of the alleles uses in the call
        alleleMap[k] = uniqueAlleles.computeIfAbsent(alleles[j].substring(leftClip), k1 -> uniqueAlleles.size());
      }
    }
    // Deal with remaining description alleles
    for (int k = 0; k < oldDescription.size(); ++k) {
      if (alleleMap[k] == -1) {
        // This description allele is not used in the actual call, it should appear in the extra alleles
        final String a = oldDescription.name(k);
        // Find in original extended alleles
        final int j = findIfPresent(extraDescriptionAlleles, a);
        assert j >= 0; // i.e. we should find it
        final List<Pair<Integer, String[]>> fullPartition = partitions.get(j + 1); // The partition corresponding to this extra allele
        final List<Pair<Integer, String[]>> partition = SplitAlleles.removeAllRef(fullPartition);
        final int offset = part.getA();
        boolean found = false;
        for (final Pair<Integer, String[]> p : partition) {
          if (p.getA() == offset) {
            // Found a piece of the partition at the same offset as the main partition
            final String[] b = p.getB();
            final String remapPart = b[b.length - 1]; // Extra allele is the last column
            if (remapPart.length() >= leftClip) {
              // Look it up in the output alleles
              final Integer outId = uniqueAlleles.get(remapPart.substring(leftClip));
              if (outId != null) {
                alleleMap[k] = outId;
                found = true;
                break;
              }
            }
          }
        }
        if (!found) {
          // Retain as is without any kind of splitting
          alleleMap[k] = uniqueAlleles.computeIfAbsent(a, k1 -> uniqueAlleles.size());
        }
      }
    }
    final String[] newUniqueAlleles = uniqueAlleles.keySet().toArray(new String[uniqueAlleles.size()]);
    return new Pair<>(newUniqueAlleles, alleleMap);
  }

  private VariantSample[] createVariantSamples(final Variant original, final SplitAlleles splitter, final String[] extraDescriptionAlleles, final List<List<Pair<Integer, String[]>>> partitions, final int slice, final int leftClip) {
    final VariantSample[] newSamples = new VariantSample[original.getNumberOfSamples()];
    final String[] originalAlleles = extractAlleles(original);
    for (int k = 0; k < newSamples.length; ++k) {
      final VariantSample sample = original.getSample(k);
      if (sample == null) {
        newSamples[k] = null;
      } else {
        final Pair<Integer, String[]> part = partitions.get(0).get(slice);
        final String newGenotype = createVariantGenotype(splitter, part.getB(), sample.getName(), leftClip);
        newSamples[k] = new VariantSample(sample.getPloidy(), newGenotype, sample.isIdentity(), sample.getMeasure(), sample.isDeNovo(), sample.getDeNovoPosterior());
        VariantSample.copy(sample, newSamples[k]);
        // Update counts to correspond to the new description
        final Description oldDescription = original.getSample(k).getStats().counts().getDescription();
        final Pair<String[], int[]> alleleMap = getAlleleMap(oldDescription, originalAlleles, extraDescriptionAlleles, partitions, slice, leftClip);
        final Description d = newSamples[k].getStats().counts().getDescription() instanceof DescriptionNone ? DescriptionNone.SINGLETON : new DescriptionCommon(alleleMap.getA());
        newSamples[k].getStats().remapAlleleStatistics(d, alleleMap.getB());
        final String va = findVariantAllele(newSamples[k], part.getB()[0], mVariantAlleleTrigger);
        newSamples[k].setVariantAllele(va == null ? null : va.substring(leftClip));
        newSamples[k].setGenotypeLikelihoods(newGenotypeLikelihoods(sample, splitter, part.getB(), leftClip));
      }
    }
    return newSamples;
  }

  private void updatePossibleCause(final Variant original, final Variant result, final SplitAlleles splitter, final String[] alleles, final int leftClip) {
    final String possibleCause = original.getPossibleCause();
    result.setPossibleCause(possibleCause == null ? null : alleles[splitter.getColumnIndex(possibleCause)].substring(leftClip));
  }

  private Variant createVariant(final Variant original, final SplitAlleles splitter, final String[] extraDescriptionAlleles, final List<List<Pair<Integer, String[]>>> partitions, final int slice, final int splitId) {
    // Get all the alleles (including reference) that are in this partition.
    // Determine the length of any common prefix (this can happen if they are
    // different lengths and depending on penalties etc. during alignment).
    // Determine the new locus represented by (the possibly clipped) alleles.
    // Construct new variant samples (i.e. genotypes) for each sample and
    // update the VAF allele (possible cause).  Update the de novo status
    // for each sample.
    final Pair<Integer, String[]> part = partitions.get(0).get(slice);
    final String[] alleles = part.getB();
    final int leftClip = StringUtils.longestPrefix(alleles);
    final VariantLocus locus = original.getLocus();
    final int offset = part.getA() + leftClip;
    final int start = locus.getStart() + offset;
    final String newRef = alleles[SplitAlleles.REFERENCE_COLUMN_INDEX].substring(leftClip);
    final int end = start + newRef.length();
    final VariantLocus newLocus = new VariantLocus(locus.getSequenceName(), start, end, newRef, getAnchorBase(locus, offset));
    final VariantSample[] newSamples = createVariantSamples(original, splitter, extraDescriptionAlleles, partitions, slice, leftClip);
    final Variant decomposedVariant = new Variant(newLocus, newSamples);
    Variant.copy(original, decomposedVariant);
    updatePossibleCause(original, decomposedVariant, splitter, alleles, leftClip);
    decomposedVariant.setSplitId(splitId);
    return mDenovoChecker != null ? denovoCorrect(mDenovoChecker, decomposedVariant) : decomposedVariant;
  }

  private String[] getExtraDescriptionAlleles(final String ref, final Set<String> alts, final Variant original) {
    final Set<String> extra = new HashSet<>();
    for (int j = 0; j < original.getNumberOfSamples(); ++j) {
      // Collect alleles from all descriptions (usally all samples have the same description, but some might be DescriptionNone)
      final VariantSample sample = original.getSample(j);
      if (sample != null) {
        final Description desc = sample.getStats().counts().getDescription();
        for (int k = 0; k < desc.size(); ++k) {
          extra.add(desc.name(k));
        }
      }
    }
    extra.removeAll(alts);
    extra.remove(ref);
    return extra.toArray(new String[extra.size()]);
  }

  @Override
  public List<Variant> decompose(final Variant original) {
    final VariantLocus locus = original.getLocus();
    final String ref = locus.getRefNts();
    final Set<String> alts = extractAlts(original);
    final String[] extraDescriptionAlleles = getExtraDescriptionAlleles(ref, alts, original);
    final SplitAlleles splitter = new SplitAlleles(ref, alts);
    final List<List<Pair<Integer, String[]>>> partitions = SplitAlleles.removeAllRefList(splitter.partition(extraDescriptionAlleles));
    final List<Variant> result = new ArrayList<>();
    int splitId = 0;
    for (int slice = 0; slice < partitions.get(0).size(); ++slice) {
      result.add(createVariant(original, splitter, extraDescriptionAlleles, partitions, slice, splitId++));
    }
    return result;
  }
}
