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
package com.rtg.ml;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Properties;

import com.rtg.util.ContingencyTable;
import com.rtg.util.IORunnable;
import com.rtg.util.PortableRandom;
import com.rtg.util.Seedable;
import com.rtg.util.SimpleThreadPool;
import com.rtg.util.ThreadAware;
import com.rtg.util.Utils;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * An ensemble classifier.
 *
 */
public class BaggedClassifierBuilder implements BuildClassifier, ThreadAware, Seedable {

  /** All our properties start with this */
  public static final String PROP_PREFIX = "bagging.";

  /** Property to set the number of sub classifiers */
  public static final String PROP_NUMTREES = PROP_PREFIX + "num-trees";
  /** Property to set the size of sub classifier training set, as a fraction of the full set */
  public static final String PROP_SUBSET_FRACTION = PROP_PREFIX + "subset-fraction";
  /** Property to set the random number seed */
  public static final String PROP_SEED = PROP_PREFIX + "seed";
  /** Property to set the implementation of sub classifier */
  public static final String PROP_SUBCLASSIFIER = PROP_PREFIX + "subclassifier";

  /** Property to set whether out of bag accuracy should be estimated */
  public static final String PROP_EVALUATE_OOB = PROP_PREFIX + "evaluate-oob";

  /** Property to set whether attribute importance should be estimated */
  public static final String PROP_EVALUATE_IMPORTANCES = PROP_PREFIX + "evaluate-importances";

  // Properties used during building, affect the resulting classifier
  /** Number of trees in the forest */
  private int mNumTrees = 10;
  /** Size of sub-tree training datasets, as a proportion of the input dataset size */
  private double mSubsetProportion = 1.0;
  /** Random number seed */
  private int mSeed = 42;
  /** Sub classifier builder type */
  private String mBuilderType = BuilderFactory.BuilderType.RANDOM_TREE.name();

  /** If true, compute out of bag error estimates */
  private boolean mEvaluateOob = true;
  /** If true, compute attribute importance estimates */
  private boolean mEvaluateImportances = true;

  /** Because we're a meta-classifier we keep any properties for passing down to child classifiers */
  private final Properties mProperties = new Properties();

  // The resulting classifier
  private BaggedClassifier mClassifier = null;

  private int mNumThreads = 1;

  private ContingencyTable mOobEvaluation = new SimpleEvaluation();


  @Override
  public void setProperties(Properties props) {
    mProperties.clear();
    mProperties.putAll(props);
    mBuilderType = props.getProperty(PROP_SUBCLASSIFIER, BuilderFactory.BuilderType.RANDOM_TREE.name());
    if (mBuilderType.equals(BuilderFactory.BuilderType.BAGGED.name())) {
      throw new IllegalArgumentException("Bagging cannot operate on itself");
    }
    setSeed(Integer.parseInt(props.getProperty(PROP_SEED, "42")));
    mNumTrees = Integer.parseInt(props.getProperty(PROP_NUMTREES, "75"));
    if (mNumTrees < 1) {
      throw new IllegalArgumentException("Number of trees must be positive");
    }
    mSubsetProportion = Double.parseDouble(props.getProperty(PROP_SUBSET_FRACTION, String.valueOf(1.0)));
    if (mSubsetProportion <= 0 || mSubsetProportion > 1.0) {
      throw new IllegalArgumentException("Subset proportion must be between 0 and 1");
    }
    mEvaluateImportances = Boolean.parseBoolean(props.getProperty(PROP_EVALUATE_IMPORTANCES, Boolean.TRUE.toString()));
    mEvaluateOob = mEvaluateImportances || Boolean.parseBoolean(props.getProperty(PROP_EVALUATE_OOB, Boolean.TRUE.toString()));
  }

  @Override
  public void setNumberOfThreads(int n) {
    mNumThreads = n;
  }

  @Override
  public int getNumberOfThreads() {
    return mNumThreads;
  }

  @Override
  public int getSeed() {
    return mSeed;
  }

  @Override
  public void setSeed(int seed) {
    mSeed = seed;
  }

  @Override
  public void build(final Dataset dataset) {

    final int numInsts = dataset.size();
    final int subsetInst = (int) (numInsts * mSubsetProportion);
    final PortableRandom random = new PortableRandom(mSeed);
    final Attribute[] attributes = dataset.getAttributes();
    final ContingencyTable[] attEvals = new ContingencyTable[attributes.length];
    if (mEvaluateImportances) {
      for (int k = 0; k < attEvals.length; ++k) {
        attEvals[k] = new SimpleEvaluation();
      }
    }

    // Build individual models
    final SimpleThreadPool stp = new SimpleThreadPool(mNumThreads, "bagger-build", false);
    final BuildClassifier[] subbuilders = new BuildClassifier[mNumTrees];
    final PredictClassifier[] subpredictors = new PredictClassifier[mNumTrees];
    final PortableRandom[] randoms = new PortableRandom[mNumTrees];
    final ContingencyTable totalEval = new SimpleEvaluation();

    for (int i = 0; i < mNumTrees; ++i) {
      final int j = i;

      randoms[j] = new PortableRandom(random.nextLong() * (j + 1) + j); // Ugly split generator
      subbuilders[j] = BuilderFactory.create(mBuilderType, mProperties);
      if (subbuilders[j] instanceof Seedable) {
        ((Seedable) subbuilders[j]).setSeed(random.nextInt() + j);  // Override seed
      }

      stp.execute(new IORunnable() {
        @Override
        public void run() {
          // Building
          final TrainTestSplit subset = TrainTestSplit.sampleWithReplacement(dataset, subsetInst, randoms[j]);
          subbuilders[j].build(subset.mTrain);
          subpredictors[j] = subbuilders[j].getClassifier();
          // Evaluating
          final Dataset testSet = subset.mTest;
          if (mEvaluateOob) {
            evaluate(totalEval, subpredictors[j], testSet); // Holdout
          }
          if (mEvaluateImportances) {
            for (int a = 0; a < attributes.length; ++a) {
              final Dataset permuted = testSet.deepCopy();
              permuteAttribute(permuted.getInstances(), a, randoms[j]);
              evaluate(attEvals[a], subpredictors[j], permuted); // Importances
            }
          }
        }
      });
    }
    try {
      stp.terminate();
    } catch (final IOException e) {
      throw new RuntimeException("Model building should not throw IOException", e);
    }
    mClassifier = new BaggedClassifier(subpredictors);
    if (mEvaluateOob) {
      mOobEvaluation = totalEval;
      final double percent = 100.0 * mOobEvaluation.accuracy();
      Diagnostic.info("Hold-out score: " + Utils.realFormat(percent, 4) + "% (" + mOobEvaluation.correct() + "/" + mOobEvaluation.total() + ")");
    }
    if (mEvaluateImportances) {
      // Calculate attribute scores
      for (int a = 0; a < attributes.length; ++a) {
        final int extraCorrect = (int) (mOobEvaluation.correct() - attEvals[a].correct());
        Diagnostic.info("Attribute importance estimate for " + attributes[a].getName() + ": " + Utils.realFormat(100.0 * extraCorrect / mOobEvaluation.correct(), 4)
            + "% (" + extraCorrect + "/" + mOobEvaluation.correct() + ")");
      }
    }
  }

  private void evaluate(ContingencyTable totalEval, PredictClassifier submodel, Dataset testSet) {
    final SimpleEvaluation eval = new SimpleEvaluation();
    eval.evaluate(submodel, testSet);
    synchronized (totalEval) {
      totalEval.add(eval);
    }
  }

  /* Randomly permute the attribute across the given set of instances */
  private static void permuteAttribute(ArrayList<Instance> instances, int attribute, PortableRandom random) {
    for (int i = instances.size() - 1; i > 0; --i) {
      final double[] current = instances.get(i).instance();
      final double[] target = instances.get(random.nextInt(i + 1)).instance();

      final double tmp = target[attribute];
      target[attribute] = current[attribute];
      current[attribute] = tmp;
    }
  }

  public double getOutOfBagAccuracy() {
    return mOobEvaluation.accuracy();
  }

  @Override
  public PredictClassifier getClassifier() {
    return mClassifier;
  }

}
