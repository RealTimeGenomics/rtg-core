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

package com.rtg.variant.bayes.multisample.population;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.launcher.globals.CoreGlobalFlags;
import com.rtg.launcher.globals.GlobalFlags;
import com.rtg.mode.DNA;
import com.rtg.tabix.TabixIndexer;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.intervals.ReferenceRanges;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.bayes.AlleleBalanceProbability;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.ModelFactory;
import com.rtg.variant.bayes.complex.DescriptionComplex;
import com.rtg.variant.bayes.multisample.HaploidDiploidHypotheses;
import com.rtg.variant.bayes.snp.DescriptionSnp;
import com.rtg.variant.bayes.snp.HypothesesNone;
import com.rtg.variant.bayes.snp.HypothesesPrior;
import com.rtg.variant.bayes.snp.HypothesesSnp;
import com.rtg.variant.bayes.snp.ModelSnpFactory;
import com.rtg.variant.util.VariantUtils;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 * Class loads a population VCF and creates hypotheses for positions in VCF file.
 */
public class PopulationHwHypothesesCreator implements SiteSpecificPriors {

  /** If true, adjust complex priors according to SSP counts, and potentially novelty. */
  static final boolean COMPLEX_SSP = true; //Boolean.valueOf(System.getProperty("rtg.ssp.complex", "true"));

  private static final int MAX_VARIANT_LENGTH = 100; //Integer.valueOf(System.getProperty("rtg.ssp.complexlength", "100"));

  // Ignore complex triggering from population variants where alternative alleles are rare
  private static final double MIN_NON_REF_FREQUENCY = GlobalFlags.getDoubleValue(CoreGlobalFlags.VARIANT_POPULATION_PRIORS_MIN_AF);

  private final HashMap<String, HashMap<Integer, HaploidDiploidHypotheses<HypothesesPrior<Description>>>> mMap = new HashMap<>();

  private final HashMap<String, ArrayList<AlleleCounts>> mAlleleCounts = new HashMap<>();

  private final Description mDescription;
  private final ModelFactory<Description, HypothesesSnp> mHaploid;
  private final ModelFactory<Description, HypothesesSnp> mDiploid;

  private final int mMaxRefLength;

  /**
   * @param input file containing population
   * @param haploidFactory factory for haploid calls
   * @param diploidFactory factory for diploid calls
   * @param ranges region restrictions
   * @exception IOException if error
   */
  PopulationHwHypothesesCreator(final File input, final ModelSnpFactory haploidFactory, final ModelSnpFactory diploidFactory, ReferenceRanges<String> ranges) throws IOException {
    mHaploid = haploidFactory;
    mDiploid = diploidFactory;
    mDescription = DescriptionSnp.SINGLETON;
    if (!TabixIndexer.indexFileName(input).exists()) {
      throw new IllegalArgumentException("Requires tabix indexed input file");
    }
    mMaxRefLength = load(input, ranges);
  }

  /**
   * @param input input file containing population
   * @param genomePriorParams genome prior params
   * @param ranges region restrictions
   * @param alleleBalance allele balance probability calculator
   * @exception IOException if error
   */
  public PopulationHwHypothesesCreator(File input, GenomePriorParams genomePriorParams, ReferenceRanges<String> ranges, AlleleBalanceProbability alleleBalance) throws IOException {
    this(input, new ModelSnpFactory(genomePriorParams, true, alleleBalance), new ModelSnpFactory(genomePriorParams, false, alleleBalance), ranges);
  }


  private int load(File inputFile, ReferenceRanges<String> ranges) throws IOException {
    int maxRefLength = 0;
    Diagnostic.developerLog("Loading population-prior allele counts...");

    try (AlleleCountsFileReader acr = AlleleCountsFileReader.openAlleleCountReader(inputFile, ranges)) {
      int totalCount = 0;
      int retainedCount = 0;
      int simpleCount = 0;
      int skippedDueToLength = 0;

      final HwEstimator hwEstimator = new HwEstimator();

      while (acr.next()) {
        final AlleleCounts alleleCounts = acr.getCurrent();
        ++totalCount;

        final String referenceName = acr.getCurrentReference();

        boolean validLength = true;
        for (String allele : alleleCounts.allelesSeen()) {
          validLength &= allele.length() <= MAX_VARIANT_LENGTH;
        }

        maxRefLength = Math.max(maxRefLength, alleleCounts.refLength());
        if (!validLength) {
          ++skippedDueToLength;
          continue;
        }

        if (COMPLEX_SSP) {
          if ((1.0 - alleleCounts.refAlleleFrequency()) >= MIN_NON_REF_FREQUENCY) {
            //System.err.println("total alleles: " + totalCount2 + ", simple? " + simple);
            final ArrayList<AlleleCounts> list = mAlleleCounts.computeIfAbsent(referenceName, k -> new ArrayList<>());
            list.add(alleleCounts);
            ++retainedCount;
          }
        }

        // simple snps from here on
        if (!alleleCounts.isComplex()) {
          final DNA reference = DNA.valueOf(alleleCounts.getReferenceAllele());
          if (reference == DNA.N) {
            continue;
          }

          final SnpDescriptionCounts c = new SnpDescriptionCounts(mDescription.size(), reference);
          for (String allele : alleleCounts.allelesSeen()) {
            final int count = alleleCounts.count(allele);
            c.increment(allele.charAt(0), count); //must be 1 length to be in here anyway.
          }

          final HaploidDiploidHypotheses<HypothesesPrior<Description>> newHyp = createHapDipHypotheses(HypothesesNone.SINGLETON, mHaploid.defaultHypotheses(c.getReferenceIndex()), mDiploid.defaultHypotheses(c.getReferenceIndex()), c, hwEstimator);

          ++simpleCount;
          final HashMap<Integer, HaploidDiploidHypotheses<HypothesesPrior<Description>>> current;
          if (mMap.containsKey(referenceName)) {
            current = mMap.get(referenceName);
          } else {
            current = new HashMap<>();
            mMap.put(referenceName, current);
          }
          final int zPos = alleleCounts.position();
          current.put(zPos, newHyp);
        }
      }
      Diagnostic.developerLog("Population-prior loaded " + totalCount + " variants");
      Diagnostic.developerLog("Population-prior ignored " + skippedDueToLength + " variants due to length > " + MAX_VARIANT_LENGTH);
      Diagnostic.developerLog("Population-prior retained " + retainedCount + " as complex triggers (non-ref frequency >= " + MIN_NON_REF_FREQUENCY + ")");
      Diagnostic.developerLog("Population-prior retained " + simpleCount + " as simple priors");
    }
    return maxRefLength;
  }

  /**
   * Create description counts based on the matches in <code>description</code>
   * @param description the description
   * @param reference the reference index
   * @return counts
   */
  public static DescriptionCounts getDescriptionCounts(DescriptionComplex description, int reference) {
    final DescriptionCounts dc = new DescriptionCounts(description.size(), reference);
    for (int i = 0; i < description.size(); ++i) {
      dc.increment(i, description.match(i).alleleCount());
    }
    return dc;
  }

  /**
   * Create a <code>HaploidDiploidHypotheses</code>, using Hardy-Weinberg estimation if description counts for the
   * population have been provided.
   * @param none the none hypotheses
   * @param haploid the haploid hypotheses
   * @param diploid the diploid hypotheses
   * @param descriptionCounts the description counts for the population
   * @param hwEstimator the Hardy-Weinberg estimator to use
   * @param <D> type of description
   * @return a <code>HaploidDiploidHypotheses</code>
   */
  public static <D extends Description> HaploidDiploidHypotheses<HypothesesPrior<D>> createHapDipHypotheses(HypothesesPrior<D> none, HypothesesPrior<D> haploid, HypothesesPrior<D> diploid, DescriptionCounts descriptionCounts, HwEstimator hwEstimator) {
    final HaploidDiploidHypotheses<HypothesesPrior<D>> tmp = new HaploidDiploidHypotheses<>(none, haploid, diploid, true, descriptionCounts);

    if (descriptionCounts != null && descriptionCounts.getTotalCount() > 0) {
      return hwEstimator.computeNewPriors(tmp, descriptionCounts.getCounts(), descriptionCounts.getTotalCount());
    } else {
      return tmp;
    }
  }

  /**
   * Adjusts hypothesis priors by applying a novelty prior factor to reference vs non-reference
   * @param hypotheses hypotheses containing input priors
   * @param novelty per base probability of encountering a novel variant
   * @param <T> type of hypotheses
   * @param <D> type of description
   * @return list containing new hypothesis, first haploid, second diploid
   */
  public static <D extends Description, T extends HypothesesPrior<D>> HaploidDiploidHypotheses<HypothesesPrior<D>> computeNoveltyPriors(HaploidDiploidHypotheses<T> hypotheses, double novelty) {

    final HypothesesPrior<D> newHaploid = rescale(hypotheses.haploid(), novelty, true);
    final HypothesesPrior<D> newDiploid = rescale(hypotheses.diploid(), novelty, false);

    return new HaploidDiploidHypotheses<>(hypotheses.none(), newHaploid, newDiploid, false, hypotheses.getDescriptionCounts());
  }

  // Rescale prior probabilities toward/away from ref such that ref probability = 1 - novelty, and all non-ref hyp sum to novelty.
  private static <D extends Description, T extends HypothesesPrior<D>> HypothesesPrior<D> rescale(final T hypotheses, final double novelty, boolean haploid) {
    final double newRefProb = 1.0 - novelty;
    final double[] priorProb = new double[hypotheses.size()];
    final PossibilityArithmetic arith = hypotheses.arithmetic();
    final int reference = hypotheses.reference();
    double total = 0;
    // First pass pull out the probabilities and total so we can normalize
    // If we know the priors are already normalized we can skip this step and maybe operate directly in possibility space
    for (int i = 0; i < priorProb.length; ++i) {
      priorProb[i] = arith.poss2Prob(hypotheses.p(i));
      total += priorProb[i];
    }
    final double hapCorrection = novelty / (total - priorProb[reference]);
    for (int i = 0; i < priorProb.length; ++i) {
      priorProb[i] = (i == reference) ? newRefProb : (hapCorrection * priorProb[i]);
    }
    return new HypothesesPrior<>(hypotheses.description(), arith, VariantUtils.prob2Poss(priorProb, arith), haploid, reference);
  }

  @Override
  public HaploidDiploidHypotheses<HypothesesPrior<Description>> getSnpHypotheses(final String templateName, int zeroPos) {
    final HashMap<Integer, HaploidDiploidHypotheses<HypothesesPrior<Description>>> hashMap = mMap.get(templateName);
    return hashMap == null ? null : hashMap.get(zeroPos);
  }

  @Override
  public List<AlleleCounts> getCounts(final String templateName, int start, int end) {
    final List<AlleleCounts> sloppyCounts = getSloppyCounts(templateName, Math.max(0, start - mMaxRefLength), end);
    sloppyCounts.removeIf(ac -> ac.position() + ac.refLength() <= start);
    return sloppyCounts;
  }

  private List<AlleleCounts> getSloppyCounts(final String templateName, int start, int end) {
    final List<AlleleCounts> acs = new ArrayList<>();

    if (mAlleleCounts.containsKey(templateName)) {
      final ArrayList<AlleleCounts> list = mAlleleCounts.get(templateName);
      int min = 0;
      int max = list.size();

      int pmin = -1;
      int pmax = -1;
      while (max > min && (min != pmin || max != pmax)) {
        final int mid = (max + min) / 2;
        pmin = min;
        pmax = max;
        final AlleleCounts ac = list.get(mid);
        final int pos = ac.position();
        if (pos < start) {
          min = mid;
        } else if (pos >= end) {
          max = mid;
        } else { // within start to end so find bounds
          min = mid - 1;
          while (min >= 0 && list.get(min).position() >= start) {
            --min;
          }
          max = mid + 1;
          while (max < list.size() && list.get(max).position() < end) {
            ++max;
          }
          for (int i = min + 1; i < max; ++i) {
            acs.add(list.get(i));
          }
          return acs;
        }
      }
    }
    return acs;
  }
}

// Keeps per base counts in 0-3 space, same as SNP hypotheses
@TestClass("com.rtg.variant.bayes.multisample.population.PopulationHwHypothesesCreatorTest")
class SnpDescriptionCounts extends DescriptionCounts {

  static final int A = DNA.A.ordinal() - 1;
  static final int C = DNA.C.ordinal() - 1;
  static final int G = DNA.G.ordinal() - 1;
  static final int T = DNA.T.ordinal() - 1;

  SnpDescriptionCounts(final int size, final DNA reference) {
    super(size, reference.ordinal() - 1); //change from 1- 4 to 0 - 3 space
  }

  void increment(char allele, int count) {
    final int index;
    switch (allele) {
      case 'a':
      case 'A' :
        index = A;
        break;

      case 'c':
      case 'C' :
        index = C;
        break;

      case 'g':
      case 'G' :
        index = G;
        break;

      case 't':
      case 'T' :
        index = T;
        break;

      default:
        index = -1;
        break;
    }
    if (index == -1) {
      return;
    }
    increment(index, count);
  }

  @Override
  public String toString() {
    return "A = " + getCount(A) + " C = " + getCount(C) + " G = " + getCount(G) + " T = " + getCount(T);
  }
}

