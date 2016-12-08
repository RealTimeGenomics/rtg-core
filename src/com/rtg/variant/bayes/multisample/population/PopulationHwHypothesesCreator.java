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

package com.rtg.variant.bayes.multisample.population;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import com.reeltwo.jumble.annotations.TestClass;
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
    if (!new File(input.getPath() + TabixIndexer.TABIX_EXTENSION).exists()) {
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
      int loadedCount = 0;
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
          //System.err.println("total alleles: " + totalCount2 + ", simple? " + simple);
          ArrayList<AlleleCounts> list = mAlleleCounts.get(referenceName);
          if (list == null) {
            list = new ArrayList<>();
            mAlleleCounts.put(referenceName, list);
          }
          list.add(alleleCounts);
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

          ++loadedCount;
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
      Diagnostic.developerLog("Population-prior loading skipped " + skippedDueToLength + " variants due to length > " + MAX_VARIANT_LENGTH);
      Diagnostic.developerLog("Population-prior loaded " + loadedCount + " / " + totalCount + " variants");
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
    for (final Iterator<AlleleCounts> it = sloppyCounts.iterator(); it.hasNext();) {
      final AlleleCounts ac = it.next();
      if (ac.position() + ac.refLength() <= start) {
        it.remove();
      }
    }
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

