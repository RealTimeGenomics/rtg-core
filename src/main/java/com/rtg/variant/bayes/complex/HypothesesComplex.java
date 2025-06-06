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

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.stream.Collectors;

import com.rtg.launcher.globals.CoreGlobalFlags;
import com.rtg.launcher.globals.GlobalFlags;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.VariantParams;
import com.rtg.variant.bayes.Code;
import com.rtg.variant.bayes.CodeDiploid;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.Model;
import com.rtg.variant.bayes.multisample.population.SiteSpecificPriors;
import com.rtg.variant.bayes.snp.HypothesesPrior;
import com.rtg.variant.match.AlignmentMatch;
import com.rtg.variant.match.AlleleAsReadMatch;
import com.rtg.variant.match.Match;
import com.rtg.variant.match.OrdinalMatchComparator;
import com.rtg.variant.util.VariantUtils;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 * Contains hypothesis for complex caller
 */
public class HypothesesComplex extends HypothesesPrior<DescriptionComplex> {

  /** print complex hypotheses for debugging */
  private static final boolean PRINT_HYP_DETAILS = GlobalFlags.isSet(CoreGlobalFlags.COMPLEX_HYPOTHESIS_DETAILS);

  private static final boolean NEW_DIPLOID_PRIORS = GlobalFlags.getBooleanValue(CoreGlobalFlags.COMPLEX_HYPOTHESIS_NEW_PRIORS);

  /** Use the genome priors and the following two constants to adjust the priors */
  private static final boolean ADJUST = GlobalFlags.getBooleanValue(CoreGlobalFlags.COMPLEX_HYPOTHESIS_ADJUST_PRIORS);

  /** Bias hypothesis priors between ref and alt alleles. 1.0 = no bias. 0.0 = reduced alt likelihood */
  private static final double PRIORS_ALT_BIAS = 0.1; //Double.parseDouble(System.getProperty("rtg.hypoth-cx-alt-bias", "0.1"));

  /** Bias hypothesis priors between heterozygous vs homozygous. 1.0 = no bias. 0.0 = reduced het likelihood */
  private static final double PRIORS_HET_BIAS = 0.01; //Double.parseDouble(System.getProperty("rtg.hypoth-cx-het-bias", "0.01"));

  /** Threshold above which hypothesis pruning will occur. Pruning will attempt to keep this many hypotheses. */
  private static final int MIN_HYPOTH = 6; //Integer.parseInt(System.getProperty("rtg.min-hypoth", "6"));
  private static final int HYPOTH_CUTOFF_DIV = 6; //Integer.parseInt(System.getProperty("rtg.hypoth-cutoff-div", "6"));
  private static final double HYPOTH_COUNTS_MULT = 1.5; //Double.parseDouble(System.getProperty("rtg.hypoth-counts-mult", "1.5"));

  /**
   * Obtain hypotheses for a complex region
   * @param context contains complex call contextual information
   * @param haploid true if reference is haploid
   * @param params variant params used for genome priors
   * @return the HypothesesComplex
   */
  public static HypothesesComplex makeComplexHypotheses(ComplexTemplate context, boolean haploid, final VariantParams params) {
    final double[] priors = makePriorsAllPaths(context, haploid, params.genomePriors());
    //assert sum >= 0 && sum <= 1;

    final HypothesesComplex result = new HypothesesComplex(context.description(), context.arithmetic(), priors, haploid, context.refHyp());

    if (PRINT_HYP_DETAILS) {
      for (int i = 0; i < result.size(); ++i) {
        final String ploidy = haploid ? "haploid" : "diploid";
        Diagnostic.developerLog("CX_PRIOR\t" + context.getSequenceName() + "\t" + context.getStart() + "\t" + context.getEnd() + "\t" + ploidy + "\t'" + result.name(i) + "'\t" + context.arithmetic().poss2Prob(result.p(i)));
      }
    }

    return result;
  }

  /**
   * Constructor that clones another HypothesesComplex but sets priors to the given priors.
   * @param hypoth HypothesesComplex to be cloned
   * @param priors new priors to use (possibility space)
   */
  public HypothesesComplex(Hypotheses<DescriptionComplex> hypoth, double[] priors) {
    super(hypoth.description(), hypoth.arithmetic(), priors, hypoth.haploid(), hypoth.reference());
  }

  HypothesesComplex(DescriptionComplex description, PossibilityArithmetic arithmetic, final double[] priors, boolean haploid, int refHyp) {
    super(description, arithmetic, priors, haploid, refHyp);
  }

  /**
   * This diagram shows the relationship of the different variables when doing the all-paths matching.
   * <img src="doc-files/makeInitialPriors.jpg" alt="image">
   * @param c region on reference being replaced by the hypotheses.
   * @param haploid true for a haploid situation
   * @param genomePriors genome prior information
   * @return the prior possibilities for transitions between haploid hypotheses.
   */
  static double[] makePriorsAllPaths(final ComplexTemplate c, boolean haploid, GenomePriorParams genomePriors) {
    if (haploid) {
      return computeHaploidPriors(c);
    } else if (NEW_DIPLOID_PRIORS) {
      return computeDiploidPriorsNew(c, genomePriors);
    } else {
      return computeDiploidPriorsOld(c, genomePriors);
    }
  }

  private static double[] computeDiploidPriorsNew(ComplexTemplate c, GenomePriorParams genomePriors) {
    final double[] haploidPriors = computeHaploidPriors(c);
    final PossibilityArithmetic arithmetic = c.arithmetic();
    final double[][] transitionProbsLn = c.transitionProbsLn();
    final int refAllele = c.refTransitionIndex();
    final Code code =  new CodeDiploid(c.description().size());
    final double[] diploidPriors = new double[code.size()];
    for (int i = 0; i < diploidPriors.length; ++i) {
      final int ca = code.a(i);
      final int cb = code.b(i);
      if (code.homozygous(i)) {
        diploidPriors[i] = haploidPriors[ca]; // prior possibility of the one allele
      } else {
        // select the most likely path to both alleles in the heterozygous case
        double p1 = arithmetic.multiply(arithmetic.ln2Poss(transitionProbsLn[refAllele][ca]), arithmetic.ln2Poss(transitionProbsLn[ca][cb]));
        final double p2 = arithmetic.multiply(arithmetic.ln2Poss(transitionProbsLn[refAllele][cb]), arithmetic.ln2Poss(transitionProbsLn[cb][ca]));
        final double p3 = arithmetic.multiply(arithmetic.ln2Poss(transitionProbsLn[refAllele][cb]), arithmetic.ln2Poss(transitionProbsLn[refAllele][ca]));
        p1 = arithmetic.gt(p1, p2) ? p1 : p2;
        p1 = arithmetic.gt(p1, p3) ? p1 : p3;
        diploidPriors[i] = arithmetic.add(p1, p1); // het genotypes require *2 factor since a|b and b|a share the same diploid Code.
      }
    }

    if (ADJUST) {
      adjustDiploidPriors(c, arithmetic, genomePriors, haploidPriors, diploidPriors);
    }
    return VariantUtils.normalisePossibilities(diploidPriors, arithmetic);
  }

  private static double[] computeDiploidPriorsOld(ComplexTemplate c, GenomePriorParams genomePriors) {
    final double[] haploidPriors = computeHaploidPriors(c);
    final PossibilityArithmetic arithmetic = c.arithmetic();
    final Code code =  new CodeDiploid(c.description().size());
    final double[] diploidPriors = new double[code.size()];
    for (int i = 0; i < diploidPriors.length; ++i) {
      final int ca = code.a(i);
      final int cb = code.b(i);
      if (code.homozygous(i)) {
        diploidPriors[i] = haploidPriors[ca]; // prior possibility of the one allele
      } else {
        double p = arithmetic.multiply(haploidPriors[ca], haploidPriors[cb]); // prior possibility of both alleles
        p = arithmetic.add(p, p); // het genotypes require *2 factor since a|b and b|a share the same diploid Code.
        diploidPriors[i] = p;
      }
    }

    if (ADJUST) {
      adjustDiploidPriors(c, arithmetic, genomePriors, haploidPriors, diploidPriors);
    }
    return VariantUtils.normalisePossibilities(diploidPriors, arithmetic);
  }

  // Adjust priors using overall genome prior information (shouldn't actually be needed as this is incorporated into all-paths alignment probabilities
  private static void adjustDiploidPriors(ComplexTemplate c, PossibilityArithmetic arithmetic, GenomePriorParams genomePriors, double[] haploidPriors, double[] diploidPriors) {
    // Reference allele frequency (average) from our priors file
    final double refFreqInitial = (1.0 - genomePriors.genomeIndelEventFraction()) / (genomePriors.genomeIndelEventFraction() + 1.0);
    double altFreqInitial = 1.0 - refFreqInitial; // Probability mass for all alt alleles.
    // TODO Hack to make alt stuff less likely, to improve ROC from slope analysis.
    altFreqInitial = altFreqInitial * PRIORS_ALT_BIAS;
    final double refFreq = 1.0 - altFreqInitial;
    final double altFreq = haploidPriors.length == 1 ? 1.0 : altFreqInitial / (haploidPriors.length - 1); // Distribute evenly among all alt alleles
    final double[] haploidFrequencies = new double[haploidPriors.length];
    for (int i = 0; i < haploidFrequencies.length; ++i) {
      haploidFrequencies[i] = arithmetic.prob2Poss(i == c.refHyp() ? refFreq : altFreq);
    }
    final double hetBias = arithmetic.prob2Poss(PRIORS_HET_BIAS);
    final Code code = new CodeDiploid(c.description().size());
    for (int i = 0; i < diploidPriors.length; ++i) {
      final int ca = code.a(i);
      final int cb = code.b(i);
      diploidPriors[i] = arithmetic.multiply(diploidPriors[i], arithmetic.multiply(haploidFrequencies[ca], haploidFrequencies[cb]));  // Probability of getting the combination of alleles (just from allele frequencies)
      if (!code.homozygous(i) && PRIORS_HET_BIAS != 1.0) {
        diploidPriors[i] = arithmetic.multiply(diploidPriors[i], hetBias); // TODO Hack to give less confidence in het stuff, to improve ROC from slope analysis.
      }
    }
  }

  // Compute prior possibilities for each description from the reference allele
  private static double[] computeHaploidPriors(ComplexTemplate c) {
    final double[] refTransitionPriorsLn = c.transitionProbsLn()[c.refTransitionIndex()];
    final PossibilityArithmetic arithmetic = c.arithmetic();
    double[] haploidPriors = new double[c.description().size()];
    for (int i = 0; i < haploidPriors.length; ++i) {
      haploidPriors[i] = arithmetic.ln2Poss(refTransitionPriorsLn[i]);
    }
    haploidPriors = VariantUtils.normalisePossibilities(haploidPriors, arithmetic);
    return haploidPriors;
  }

  /**
   * Create the complex description from the alignment matches, site-specific priors etc
   * @param readMatches all alignment matches
   * @param reference used to create reference description if not present in the matches
   * @param ssp site-specific priors for additional description entries
   * @param pruneMatches if set, prune full description
   * @param maxHypotheses maximum number of descriptions
   * @return the description
   */
  public static DescriptionComplex createComplexDescription(List<AlignmentMatch> readMatches, ComplexTemplate reference, SiteSpecificPriors ssp, boolean pruneMatches, int maxHypotheses) {

    // turn ssps into alleleAsRef matches
    final ArrayList<Match> sspMatches = new ArrayList<>();
    // TODO uncomment the lines below once we figured out how to create alleles based on phasing information
    //      currently adding any alleles with their respective counts creates issues for low coverage data
    //      for example, A = 0, C = 1000, T = 1100, G = 0, if we try to make a call where reference in T and
    //      all the low coverage evidence suggest C, the correct call would have been a C, but adding
    //      those alleles will make the call T:C
    //    int refCount = 0;
    //    if (ssp != null) {
    //      // add allele count matches
    //      final List<AlleleCounts> counts = ssp.getCounts(reference.refName(), reference.start(), reference.end());
    //      if (counts.size() > 1) {
    //        Diagnostic.developerLog("loading " + counts.size() + " complex SSPs for " + reference.refName() + ":" + reference.start() + "-" + reference.end());
    //      }
    //      for (final AlleleCounts ac : counts) {
    //        for (final String allele : ac.allelesSeen()) {
    //          final byte[] alleleBytes = refSeqFromAllele(reference, ac.position(), ac.refLength(), allele);
    //          if (!Arrays.equals(alleleBytes, reference.replaceBytes())) {
    //            // don't want to add unchanged reference to set of matches
    //            sspMatches.add(new AlleleAsReadMatch(alleleBytes, ac.count(allele)));
    //          } else {
    //            refCount += ac.count(allele);
    //          }
    //        }
    //      }
    //    }
    //    final Match ref = new AlleleAsReadMatch(reference.replaceBytes(), refCount);

    // Extract unique hypotheses from the matches, discarding those that contain
    // an N or do not span the entire region of interest.
    final List<AlignmentMatch> filteredMatches = readMatches.stream().filter(m -> m.mapError() < Model.AMBIGUITY_THRESHOLD && !m.readString().contains("N") && m.isFixedLeft() && m.isFixedRight()).collect(Collectors.toList());
    final Match ref = new AlleleAsReadMatch(reference.replaceBytes());
    final ArrayList<Match> res;
    if (pruneMatches) { //TODO eventually, do this in all cases (if deemed acceptable)
      res = createPrunedDescription(filteredMatches, ref, sspMatches, maxHypotheses);
    } else {
      // TreeSet is used to ensure a consistent ordering of hypotheses in the result.
      final TreeSet<Match> matchSet = new TreeSet<>(new OrdinalMatchComparator());
      if (!ref.readString().contains("N")) {
        matchSet.add(ref);
      }
      matchSet.addAll(filteredMatches);
      res = new ArrayList<>(matchSet);
      res.addAll(sspMatches);
    }
    return new DescriptionComplex(res);
  }

  /**
   * Sort by frequency of match descending then in {@code OrdinalMatchComparator}
   */
  private static class MatchCountComparator implements Comparator<MatchCount>, Serializable {
    private static final OrdinalMatchComparator COMPARE = new OrdinalMatchComparator();

    @Override
    public int compare(MatchCount o1, MatchCount o2) {
      if (o1 == null) {
        if (o2 == null) {
          return 0;
        }
        return -1;
      } else if (o2 == null) {
        return 1;
      }
      final int count = Integer.compare(o2.count(), o1.count());
      if (count != 0) {
        return count;
      }
      return COMPARE.compare(o1.match(), o2.match());


    }
  }

  private static class MatchCount {
    final Match mMatch;
    private int mCount;
    private int mId = -1;
    private int mVaf = 0;
    MatchCount(Match match) {
      mMatch = match;
    }
    Match match() {
      return mMatch;
    }
    void increment() {
      ++mCount;
    }
    int count() {
      return mCount;
    }
    void setId(int id) {
      mId = id;
    }
    int getId() {
      return mId;
    }
    void setVaf(int vaf) {
      mVaf = vaf;
    }
    int getVaf() {
      return mVaf;
    }
  }

  private static ArrayList<Match> createPrunedDescription(List<AlignmentMatch> matches, Match ref, ArrayList<Match> sspMatches, int maxHypotheses) {
    // TreeSet is used to ensure a consistent ordering of hypotheses in the result.
    final TreeMap<String, MatchCount> matchMap = new TreeMap<>();

    final boolean refOk = !ref.readString().contains("N");
    if (refOk) {
      matchMap.put(ref.readString(), new MatchCount(ref));
    }

    for (final Match match : sspMatches) {
      matchMap.put(match.readString(), new MatchCount(match));
    }

    // Extract unique hypotheses from the matches
    int maxGenome = 0;
    for (final AlignmentMatch m : matches) {
      maxGenome = Math.max(maxGenome, m.alignmentRecord().getGenome());
      matchMap.computeIfAbsent(m.readString(), ms -> new MatchCount(m)).increment();
    }

    if (matchMap.size() == 0) {
      return new ArrayList<>();
    }

    final List<MatchCount> sorted = new ArrayList<>(matchMap.values());
    sorted.sort(new MatchCountComparator());

    final int size = sorted.size();
    final int highCount = sorted.get(0).count();
    final int lowCount = sorted.get(size - 1).count();
    if (matchMap.size() > 0) {
      //filter low count hypotheses if there are sufficient different hypotheses
      //use the last hypothesis to determine cutoff value ( div 6? )
      //keep at least 6 hypotheses then find the cutoff value below that (taking into account the value of the 6th hypothesis)
      //if the 6th hypothesis is below the cutoff pick a cutoff that will include it.
      final int cutoff;
      if (size > MIN_HYPOTH) {
        final int cutoff1 = lowCount / HYPOTH_CUTOFF_DIV;
        final int totalCoverage = matches.size();
        if (sorted.get(size - MIN_HYPOTH).count() == 1 && totalCoverage > size * HYPOTH_COUNTS_MULT) {
          cutoff = 2;
        } else if (sorted.get(size - MIN_HYPOTH).count() < cutoff1) {
          cutoff = sorted.get(size - MIN_HYPOTH).count();
        } else {
          cutoff = cutoff1;
        }
      } else {
        cutoff = Integer.MIN_VALUE;
      }
      //System.err.println("cutoff: " + cutoff);
//      for (MatchCount mc : sorted) {
//        Diagnostic.developerLog("HypothesesComplex: " + (mc.count() < cutoff ? "removed" : "keeping") + " " + mc.match() + " count=" + mc.count());
//      }
      sorted.removeIf(m -> m.count() < cutoff);
    }

    if (sorted.size() > maxHypotheses) {
      // Sort by mean VAF in the genomes observed in, and keep the top maxHypotheses
      final int numHyps = sorted.size();
      final int numGenomes = maxGenome + 1;

      int i = 0;
      for (MatchCount mc : sorted) {
        mc.setId(i++);
      }
      final int[][] ads = new int[numGenomes][numHyps];
      final int[] dps = new int[numGenomes];
      for (final AlignmentMatch m : matches) {
        final MatchCount mc = matchMap.get(m.readString());
        if (mc.getId() != -1) {
          final int g = m.alignmentRecord().getGenome();
          ++ads[g][mc.getId()];
          ++dps[g];
        }
      }

      final int[] seenIn = new int[sorted.size()];
      final int[] sumVaf = new int[sorted.size()];
      for (int h = 0; h < numHyps; ++h) {
        for (int g = 0; g < numGenomes; ++g) {
          if (ads[g][h] > 1) {
            ++seenIn[h];
            sumVaf[h] += ads[g][h] * 100 / (dps[g] + 1);
          }
        }
        if (seenIn[h] > 0) {
          sorted.get(h).setVaf(sumVaf[h] / seenIn[h]);
        }
      }
      sorted.sort(Comparator.comparingInt(o -> -o.getVaf()));

//      int j = 0;
//      for (MatchCount mc : sorted) {
//        Diagnostic.developerLog("HypothesesComplex: " + (++j <= maxHypotheses ? "keeping" : "removed") + " " + mc.match() + " vaf=" + mc.getVaf());
//      }
    }

    final List<MatchCount> keep = sorted.subList(0, Math.min(maxHypotheses, sorted.size()));
    if (keep.size() < size) {
      Diagnostic.developerLog("Hypotheses extracted from " + matches.size() + " matches, with counts ranging from "
                              + highCount + " to " + lowCount
                              + " pruned from " + size + " to " + keep.size());
    }
    final ArrayList<Match> result = new ArrayList<>();
    for (MatchCount val : keep) {
      result.add(val.match());
    }
    if (!result.contains(ref) && refOk) {
      result.add(ref);
    }
    for (final Match match : sspMatches) { // add ssp matches if they were pruned out above
      if (!result.contains(match)) {
        result.add(match);
      }
    }
    //    System.err.println("tot " + matchMap.size());
    //    System.err.println();
    //    System.err.println();
    return result;
  }

}
