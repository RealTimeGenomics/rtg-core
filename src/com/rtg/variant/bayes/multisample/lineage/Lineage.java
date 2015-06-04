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
package com.rtg.variant.bayes.multisample.lineage;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.rtg.relation.LineageLookup;
import com.rtg.util.Pair;
import com.rtg.variant.VariantOutputLevel;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantSample;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.GenotypeMeasure;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.multisample.AbstractMultisampleCaller;
import com.rtg.variant.bayes.multisample.HaploidDiploidHypotheses;
import com.rtg.variant.bayes.snp.HypothesesPrior;

/**
 * Tree structure for a lineage.
 *
 */
public final class Lineage extends AbstractMultisampleCaller {

  /**
   * Constructs lineages
   */
  public static class LineageBuilder {

    private Map<Integer, Set<Integer>> mLineageGraph = new HashMap<>();
    private Map<Integer, Double> mDeNovoPriors = new HashMap<>();
    private double mDenovoPriorDefault = 0;
    private int mMaxId = 0;
    boolean mCoverage = false;
    VariantParams mParams;

    /**
     * Add a parent child relationship
     * @param parent the parent id
     * @param child the child id
     * @return this so calls can be chained
     */
    public LineageBuilder add(final int parent, final int child) {
      // todo should this be in terms of sample names?
      final Set<Integer> children = mLineageGraph.get(parent);
      if (children != null) {
        children.add(child);
      } else {
        final HashSet<Integer> c = new HashSet<>();
        c.add(child);
        mLineageGraph.put(parent, c);
      }
      mMaxId = Math.max(mMaxId, Math.max(parent, child));
      return this;
    }
    /**
     * Should the caller consider copy number variation
     * @param coverage true if the caller should consider copy number
     * @return this so calls can be chained
     */
    public LineageBuilder coverage(boolean coverage) {
      mCoverage = coverage;
      return this;
    }

    /**
     * The default probability of encountering a de novo mutation
     * @param prior the de novo mutation chance in raw probability space
     * @return this so calls can be chained
     */
    public LineageBuilder deNovoPriorDefault(final double prior) {
      if (prior < 0 || prior > 1) {
        throw new IllegalArgumentException();
      }
      mDenovoPriorDefault = prior;
      return this;
    }

    /**
     * Set the de novo probability for an individual sample
     * @param sample the sample id
     * @param prior the de novo mutation chance in raw probability space
     * @return this so calls can be chained
     */
    public LineageBuilder deNovoPrior(final int sample, final double prior) {
      if (prior < 0 || prior > 1) {
        throw new IllegalArgumentException();
      }
      mMaxId = Math.max(mMaxId, sample);
      mDeNovoPriors.put(sample, prior);
      return this;
    }

    /**
     * Variant caller parameters
     * @param params a params with the output options needed
     * @return this so calls can be chained
     */
    public LineageBuilder params(VariantParams params) {
      mParams = params;
      return this;
    }

    /**
     * Construct the lineage
     * @return a new lineage object as configured by this builder
     */
    public Lineage create() {
      final double[] dnp = new double[mMaxId + 1];
      Arrays.fill(dnp, mDenovoPriorDefault);
      for (Map.Entry<Integer, Double> e : mDeNovoPriors.entrySet()) {
        dnp[e.getKey()] = e.getValue();
      }
      return new Lineage(mMaxId, mLineageGraph, dnp, mCoverage, mParams);
    }
  }

  private static final int MAX_COPY_NUMBER = 2;

  private final Map<Integer, Set<Integer>> mLineageGraph;
  //private final Variable[] mGenotypeVariables;
  private final Variable[] mCoverageVariables;
  private final LineageLookup mLineageLookup;
  private final double[] mDenovoPriors;
  private final VariantParams mParams;

  private Lineage(final int max, final Map<Integer, Set<Integer>> lineageGraph, final double[] denovoPriors, boolean coverage, VariantParams params) {
    mParams = params;
    mLineageGraph = lineageGraph;
    mDenovoPriors = denovoPriors;
    mCoverageVariables = coverage ? new Variable[max + 1] : null;
    if (coverage) {
      for (int k = 0; k < mCoverageVariables.length; k++) {
        mCoverageVariables[k] = new Variable("C" + k, MAX_COPY_NUMBER + 1);
      }
    }
    final int[] parents = new int[max + 1];
    Arrays.fill(parents, -1);
    for (final Map.Entry<Integer, Set<Integer>> e : mLineageGraph.entrySet()) {
      final int parent = e.getKey();
      for (final int child : e.getValue()) {
        if (parents[child] != -1) {
          throw new IllegalArgumentException("Child " + child + " has multiple parents");
        }
        parents[child] = parent;
      }
    }
    mLineageLookup = new LineageLookup(parents);
  }

  Set<Integer> children(final int node) {
    final Set<Integer> children = mLineageGraph.get(node);
    return children == null ? Collections.<Integer>emptySet() : children;
  }

  LineageLookup lineageLookup() {
    return mLineageLookup;
  }

  boolean isRoot(final int node) {
    return mLineageLookup.getOriginal(node) == -1;
  }

  /**
   * Return parent id of a given node, or -1 if the node has no parents.
   * @param node the child node
   * @return the parent node or -1
   */
  public int parent(final int node) {
    return mLineageLookup.getOriginal(node);
  }

  Variable getCoverageVariable(final int node) {
    // todo return null if coverage not active
    return mCoverageVariables == null ? null : mCoverageVariables[node];
  }

  double deNovoPrior(final int node) {
    return mDenovoPriors[node];
  }

  @Override
  protected <D extends Description, T extends HypothesesPrior<D>> ComparisonResult makeSamples(List<ModelInterface<?>> models, HaploidDiploidHypotheses<T> hypotheses) {
    // todo Coverage
    final Factor[] rootFactors = new Factor[models.size()];
    for (int i = 0; i < models.size(); i++) {
      if (isRoot(i)) {
        final T h = hypotheses.get(models.get(i));
        rootFactors[i] = new ModelFactor(new Variable("G" + i, h.size()), h);
      }

    }
    final ForwardBackwardLineage fb = new ForwardBackwardLineage(models.get(0).arithmetic(), this, rootFactors, models);
    final VariantSample[] samples = new VariantSample[models.size()];
    boolean interesting = false;
    for (int k = 0; k < models.size(); k++) {
      final Double deNovoScore;
      final VariantSample.DeNovoStatus deNovo;
      if (isRoot(k)) {
        deNovo = VariantSample.DeNovoStatus.UNSPECIFIED;
        deNovoScore = null;
      } else {
        final Factor posteriorDeNovo = fb.posteriorDeNovo(k);
        final DefaultFactor marginal = DefaultFactor.asDefault(posteriorDeNovo.marginal(Collections.singleton(ForwardBackwardLineage.DE_NOVO)));
        final Pair<Map<Variable, Integer>, Double> bestDenovo = marginal.best();
        deNovo = bestDenovo.getA().get(ForwardBackwardLineage.DE_NOVO) == 0 ? VariantSample.DeNovoStatus.NOT_DE_NOVO : VariantSample.DeNovoStatus.IS_DE_NOVO;
        deNovoScore = marginal.arithmetic().poss2Ln(bestDenovo.getB());
      }
      final ModelInterface<?> model = models.get(k);
      final Hypotheses<? extends Description> sampleHyp = model.hypotheses();
      final DefaultFactor posterior = DefaultFactor.asDefault(fb.posterior(k));
      assert posterior.scope().size() == 1; // Should be over just genotype
      final GenotypeMeasure measure = new DefaultFactor.FactorGenotypeMeasure(posterior, model.hypotheses());
      final boolean match = measure.best() == measure.reference();

      final VariantSample sample = new VariantSample(sampleHyp.ploidy(), sampleHyp.name(measure.best()), match, measure, deNovo, deNovoScore);
      model.statistics().addCountsToSample(sample, model, mParams);
      samples[k] = sample;
      interesting |= !match;
    }
    if (!interesting && mParams.callLevel() != VariantOutputLevel.ALL) {
      return null;
    }
    return new ComparisonResult(interesting, samples, null);
  }

  @Override
  protected VariantParams getParams() {
    return mParams;
  }
}
