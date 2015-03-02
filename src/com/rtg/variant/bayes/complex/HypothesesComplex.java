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

import java.io.Serializable;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.security.AccessController;
import java.security.PrivilegedAction;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.TreeSet;

import com.rtg.mode.DNA;
import com.rtg.util.SeparateClassLoader;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.SlimException;
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
import com.rtg.variant.realign.AbstractAllPaths;
import com.rtg.variant.realign.AlignmentEnvironment;
import com.rtg.variant.realign.AlignmentEnvironmentGenome;
import com.rtg.variant.realign.AlignmentEnvironmentGenomeSubstitution;
import com.rtg.variant.realign.AllPaths;
import com.rtg.variant.realign.Environment;
import com.rtg.variant.realign.EnvironmentCombined;
import com.rtg.variant.realign.RealignParams;
import com.rtg.variant.realign.RealignParamsGenome;
import com.rtg.variant.realign.ScoreFastUnderflow;
import com.rtg.variant.realign.ScoreMatrix;
import com.rtg.variant.util.VariantUtils;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 * Contains hypothesis for complex caller
 */
public class HypothesesComplex extends HypothesesPrior<DescriptionComplex> {

  /** print complex hypotheses for debugging */
  private static final boolean PRINT_HYP_DETAILS = false; //Boolean.valueOf(System.getProperty("rtg.print-cx-hyp-details", "false"));

  /** Normalize priors after initial calculation */
  private static final boolean NORMALIZE_PRIORS = true; //Boolean.parseBoolean(System.getProperty("rtg.hypoth-cx-prior-normalisation", "true"));

  /** Bias hypothesis priors between ref and alt alleles. 1.0 = no bias. 0.0 = reduced alt likelihood */
  private static final double PRIORS_ALT_BIAS = 0.1; //Double.parseDouble(System.getProperty("rtg.hypoth-cx-alt-bias", "0.1"));

  /** Bias hypothesis priors between heterozygous vs homozygous. 1.0 = no bias. 0.0 = reduced het likelihood */
  private static final double PRIORS_HET_BIAS = 0.01; //Double.parseDouble(System.getProperty("rtg.hypoth-cx-het-bias", "0.01"));

  /** Threshold above which hypothesis pruning will occur. Pruning will attempt to keep this many hypotheses. */
  private static final int MIN_HYPOTH = 6; //Integer.parseInt(System.getProperty("rtg.min-hypoth", "6"));
  private static final int HYPOTH_CUTOFF_DIV = 6; //Integer.parseInt(System.getProperty("rtg.hypoth-cutoff-div", "6"));
  private static final double HYPOTH_COUNTS_MULT = 1.5; //Double.parseDouble(System.getProperty("rtg.hypoth-counts-mult", "1.5"));


  /*
   * The following block of code is being used to allow two separate class definitions
   * for classes that the java JIT can optimise quite heavily in another part of the code
   * if not used directly in this code. This allows it to maintain the speed it used to have
   * instead of becoming approximately 40% slower due to the optimisations removed to allow
   * both locations to use the same code.
   */
  /*******************START BLOCK******************/
  private static final Constructor<AllPaths> SCORE_FAST_CONSTRUCTOR;
  private static final Constructor<Environment> ENVIRONMENT_COMBINED_CONSTRUCTOR;
  static {
    final SeparateClassLoader loader = AccessController.doPrivileged(new PrivilegedAction<SeparateClassLoader>() {
      @Override
      public SeparateClassLoader run() {
        return new SeparateClassLoader(ScoreFastUnderflow.class, ScoreMatrix.class, AbstractAllPaths.class, EnvironmentCombined.class);
      }
    });
    try {
      @SuppressWarnings("unchecked")
      final Class<AllPaths> sfClazz = (Class<AllPaths>) loader.loadClass(ScoreFastUnderflow.class.getName());
      SCORE_FAST_CONSTRUCTOR = sfClazz.getConstructor(RealignParams.class);
      @SuppressWarnings("unchecked")
      final Class<Environment> ecClazz = (Class<Environment>) loader.loadClass(EnvironmentCombined.class.getName());
      ENVIRONMENT_COMBINED_CONSTRUCTOR = ecClazz.getConstructor(AlignmentEnvironment.class, int.class, int.class, AlignmentEnvironment.class);
    } catch (final ClassNotFoundException | NoSuchMethodException | SecurityException e) {
      throw new SlimException(e);
    }
  }

  private static AllPaths initScoreFastUnderflow(RealignParams params) {
    try {
      return SCORE_FAST_CONSTRUCTOR.newInstance(params);
    } catch (final IllegalArgumentException | InvocationTargetException | IllegalAccessException | InstantiationException e) {
      throw new SlimException(e);
    }
    //return new ScoreFastUnderflow(params);
  }

  private static Environment initEnvironmentCombined(AlignmentEnvironment samEnv, int zeroBasedStartPos, int maxShift, AlignmentEnvironment templateEnv) {
    try {
      return ENVIRONMENT_COMBINED_CONSTRUCTOR.newInstance(samEnv, zeroBasedStartPos, maxShift, templateEnv);
    } catch (final IllegalArgumentException | InvocationTargetException | IllegalAccessException | InstantiationException e) {
      throw new SlimException(e);
    }
    //return new EnvironmentCombined(samEnv, zeroBasedStartPos, maxShift, templateEnv);
  }
  /*******************END BLOCK********************/


  /**
   * Obtain hypotheses for a complex region
   * @param reference reference bases for complex region
   * @param matches matches covering complex region (see <code>intersectSet</code>)
   * @param arithmetic implementation to use for arithmetic
   * @param haploid true if reference is haploid
   * @param params variant params used for genome priors
   * @param ssp site specific priors
   * @return the HypothesesComplex
   */
  public static HypothesesComplex makeComplexHypotheses(ComplexTemplate reference, List<AlignmentMatch> matches, PossibilityArithmetic arithmetic, boolean haploid, final VariantParams params, SiteSpecificPriors ssp) {
    final DescriptionComplex description = createDescription(matches, reference, params.pruneHypotheses(), ssp, params.maxComplexHypotheses());

    int refHyp1 = -1;
    for (int i1 = 0; i1 < description.size(); i1++) {
      if (description.name(i1).equals(reference.replaceString())) {
        refHyp1 = i1;
        break;
      }
    }
    final int refHyp = refHyp1;

    double[] priors = makePriorsAllPaths(description, haploid, reference, arithmetic, refHyp, params.genomePriors());
    priors = NORMALIZE_PRIORS ? VariantUtils.normalisePossibilities(priors, arithmetic) : priors;
    //assert sum >= 0 && sum <= 1;

    final HypothesesComplex result = new HypothesesComplex(description, arithmetic, priors, haploid, refHyp);

    if (PRINT_HYP_DETAILS) {
      for (int i = 0; i < result.size(); i++) {
        final String ploidy = haploid ? "haploid" : "diploid";
        System.err.println("" + reference.getSequenceName() + ":" + reference.getStart() + "-" + reference.getEnd() + " ploidy=" + ploidy + " hyp=" + result.name(i) + " prior=" +  arithmetic.poss2Prob(result.p(i)));
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
   * @param description description for hypotheses
   * @param haploid true for a haploid situation
   * @param cot region on reference being replaced by the hypotheses.
   * @param arithmetic implementation to use for arithmetic
   * @param refHyp index of reference hypothesis
   * @param genomePriors genome prior information
   * @return the prior possibilities for transitions between haploid hypotheses.
   */
  static double[] makePriorsAllPaths(DescriptionComplex description, boolean haploid, final ComplexTemplate cot, PossibilityArithmetic arithmetic, int refHyp, GenomePriorParams genomePriors) {
    // Compute likelihood of the alleles themselves arising
    final AllPaths sm = initScoreFastUnderflow(RealignParamsGenome.SINGLETON);
    final int hypExtension = Math.max(5, Math.max(description.maxLength() + 1, cot.getLength() + 1));
    final int start = cot.getStart() - hypExtension;
    final int end = cot.getStart() + hypExtension;
    final int e0Hx = cot.getEnd() - hypExtension;
    final int refLen = cot.getLength();
    final AlignmentEnvironment[] hypEnvs = new AlignmentEnvironment[description.size()];
    final int[] hypLengths = new int[description.size()];
    final List<AlignmentEnvironment> alleles = new ArrayList<>();
    for (int i = 0; i < description.size(); i++) {
      final AlignmentEnvironment env;
      if (i == refHyp) {
        env = new AlignmentEnvironmentGenome(start, end, cot.templateBytes());
        hypLengths[i] = cot.getEnd() - cot.getStart();
      } else {
        final byte[] hypDna = DNA.stringDNAtoByte(description.name(i));
        hypLengths[i] = hypDna.length;
        env = new AlignmentEnvironmentGenomeSubstitution(start, end, cot, hypDna);
      }
      hypEnvs[i] = env;
      if (i == refHyp) {
        alleles.add(env);
      }
    }
    if (alleles.size() == 0) { // Reference wasn't included in the set of hypotheses (perhaps due to Ns)
      alleles.add(new AlignmentEnvironmentGenome(start, end, cot.templateBytes()));
    }
    double[] haploidPriors = new double[description.size()];
    Arrays.fill(haploidPriors, arithmetic.zero());
    for (final AlignmentEnvironment refEnv : alleles) {
      for (int i = 0; i < description.size(); i++) {
        final AlignmentEnvironment hypEnv = hypEnvs[i];
        final int alignStart = e0Hx - (hypLengths[i] + refLen) / 2;
        final int maxShift = (hypLengths[i] + refLen + 1) / 2;
        final Environment env = initEnvironmentCombined(hypEnv, alignStart, maxShift, refEnv);
        sm.setEnv(env);
        haploidPriors[i] = arithmetic.add(haploidPriors[i], arithmetic.ln2Poss(sm.totalScoreLn()));
      }
    }

    if (haploid) {
      // We could apply an allele frequency correction here too, but don't have this in our priors.
      return haploidPriors;
    } else {

      haploidPriors = VariantUtils.normalisePossibilities(haploidPriors, arithmetic);

      // Compute diploid priors from the haploid distances

      // Reference allele frequency (average) from our priors file
      final double refFreqInitial = (1.0 - genomePriors.genomeIndelEventFraction()) / (genomePriors.genomeIndelEventFraction() + 1.0);
      double altFreqInitial = 1.0 - refFreqInitial; // Probability mass for all alt alleles.

      // XXX Hack to make alt stuff less likely, to improve ROC from slope analysis.
      altFreqInitial = altFreqInitial * PRIORS_ALT_BIAS;

      final double refFreq = 1.0 - altFreqInitial;
      final double altFreq = haploidPriors.length == 1 ? 1.0 : altFreqInitial / (haploidPriors.length - 1); // Distribute evenly among all alt alleles

      final double[] haploidFrequencies = new double[description.size()];
      for (int i = 0; i < haploidFrequencies.length; i++) {
        haploidFrequencies[i] = arithmetic.prob2Poss(i == refHyp ? refFreq : altFreq);
      }

      final Code code =  new CodeDiploid(description.size());
      final double hetBias = arithmetic.prob2Poss(PRIORS_HET_BIAS);
      final double[] diploidPriors = new double[code.size()];
      for (int i = 0; i < diploidPriors.length; i++) {
        final int ca = code.a(i);
        final int cb = code.b(i);
        double p = arithmetic.multiply(haploidFrequencies[ca], haploidFrequencies[cb]);  // Probability of getting the combination of alleles (just from allele frequencies)
        p = arithmetic.multiply(p, haploidPriors[ca]); // And factor in the prior probability of the first allele
        if (!code.homozygous(i)) {
          p = arithmetic.multiply(p, haploidPriors[cb]); // Need to factor in the prior probability of second allele (since it is different to the first)
          p = arithmetic.add(p, p); // het genotypes have *2 factor since a|b and b|a share the same diploid Code.

          if (PRIORS_HET_BIAS != 1.0) {
            p = arithmetic.multiply(p, hetBias); // XXX Hack to give less confidence in het stuff, to improve ROC from slope analysis.
          }

        }
        diploidPriors[i] = p;
      }
      return diploidPriors;
    }
  }

  static DescriptionComplex createDescription(List<AlignmentMatch> matches, ComplexTemplate reference, boolean pruneMatches, SiteSpecificPriors ssp, int maxHypotheses) {

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

    final Match ref = new AlleleAsReadMatch(reference.replaceBytes());
    final ArrayList<Match> res;
    if (pruneMatches) { //TODO eventually, do this in all cases (if deemed acceptable)
      res = createPrunedDescription(matches, ref, sspMatches, maxHypotheses);
    } else {
      // Extract unique hypotheses from the matches, discarding those that contain
      // an N or do not span the entire region of interest. TreeSet is used to
      // ensure a consistent ordering of hypotheses in the result.
      final TreeSet<Match> matchSet = new TreeSet<>(new OrdinalMatchComparator());
      if (!ref.readString().contains("N")) {
        matchSet.add(ref);
      }
      for (final Match m : matches) {
        if (m.mapError() < Model.AMBIGUITY_THRESHOLD && !m.readString().contains("N") && m.isFixedLeft() && m.isFixedRight()) {
          matchSet.add(m);
        }
      }
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
    final int mCount;
    MatchCount(Match match, int count) {
      mMatch = match;
      mCount = count;
    }
    Match match() {
      return mMatch;
    }

    int count() {
      return mCount;
    }

  }

  private static ArrayList<Match> createPrunedDescription(List<AlignmentMatch> matches, Match ref, ArrayList<Match> sspMatches, int maxHypotheses) {
    final TreeMap<Match, Integer> matchMap = new TreeMap<>(new OrdinalMatchComparator());

    final boolean refOk = !ref.readString().contains("N");
    if (refOk) {
      matchMap.put(ref, 0);
    }

    for (final Match match : sspMatches) {
      matchMap.put(match, 0);
    }

    // Extract unique hypotheses from the matches, discarding those that contain
    // an N or do not span the entire region of interest. TreeSet is used to
    // ensure a consistent ordering of hypotheses in the result.
    int totalCoverage = 0;
    for (final Match m : matches) {
      if (m.mapError() < Model.AMBIGUITY_THRESHOLD && !m.readString().contains("N") && m.isFixedLeft() && m.isFixedRight()) {
        if (matchMap.containsKey(m)) {
          int i = matchMap.get(m);
          matchMap.put(m, ++i);
        } else {
          matchMap.put(m, 1);
        }
        totalCoverage++;
      }
    }

    if (matchMap.size() == 0) {
      return new ArrayList<>();
    }

    final List<MatchCount> sorted = new ArrayList<>();
    for (Entry<Match, Integer> entry : matchMap.entrySet()) {
      sorted.add(new MatchCount(entry.getKey(), entry.getValue()));
    }
    Collections.sort(sorted, new MatchCountComparator());

    final int size = sorted.size();
    final int biggestVal = sorted.get(size - 1).count();
    if (matchMap.size() > 0) {
      //filter low count hypotheses if there are sufficient different hypotheses
      //use the last hypothesis to determine cutoff value ( div 6? )

      int cutoff = biggestVal / HYPOTH_CUTOFF_DIV;

      //keep at least 6 hypotheses then find the cutoff value below that (taking into account the value of the 6th hypothesis)

      //if the 6th hypothesis is below the cutoff pick a cutoff that will include it.
      if (size > MIN_HYPOTH) {

        if (sorted.get(size - MIN_HYPOTH).count() == 1 && totalCoverage > size * HYPOTH_COUNTS_MULT) {
          cutoff = 2;
        } else if (sorted.get(size - MIN_HYPOTH).count() < cutoff) {
          cutoff = sorted.get(size - MIN_HYPOTH).count();
        }
      } else {
        cutoff = Integer.MIN_VALUE;
      }
      //System.err.println("cutoff: " + cutoff);

      final Iterator<MatchCount> it = sorted.iterator();
      while (it.hasNext()) {
        final MatchCount e = it.next();
        if (e.count() < cutoff) {
          //          System.err.println("removing; " + e.getKey() + " LL " + e.getValue());
          it.remove();
        } else {
          //          System.err.println("keeping: " + e.getKey() + " :: " + e.getValue());
        }
      }
    }
    final List<MatchCount> keep = sorted.subList(0, Math.min(maxHypotheses, sorted.size()));
    if (keep.size() < size) {
      Diagnostic.developerLog("Hypotheses extracted from " + matches.size() + " matches, with counts ranging from "
                              + sorted.get(0).count() + " to " + biggestVal
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
