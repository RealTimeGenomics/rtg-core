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

/**
 * Evaluate a classifier with respect to a test dataset.
 */
public class SimpleEvaluation {

  private static final int NEG = 0;
  private static final int POS = 1;

  /** First dimension is actual, second dimension is predicted */
  private final double[][] mContingencyTable = new double[2][2];

  /**
   * Evaluate the classifier on the supplied dataset
   * @param classifier the classifier to evaluate
   * @param dataset the dataset to evaluate on
   */
  public void evaluate(PredictClassifier classifier, Dataset dataset) {
    for (Instance instance : dataset.getInstances()) {
      final int actual = instance.isPositive() ? POS : NEG;
      final int predicted = classifier.predict(instance.instance()) > 0.5 ? POS : NEG;
      mContingencyTable[actual][predicted] += instance.weight();
    }
  }

  /**
   * Add in the results from the supplied evaluation
   * @param eval another evaluation
   */
  public void add(SimpleEvaluation eval) {
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        mContingencyTable[i][j] += eval.mContingencyTable[i][j];
      }
    }
  }

  /** @return the number of true positives. */
  public int truePositives() {
    return (int) mContingencyTable[POS][POS];
  }

  /** @return the number of true positives. */
  public int trueNegatives() {
    return (int) mContingencyTable[NEG][NEG];
  }

  /** @return the number of true positives. */
  public int falsePositives() {
    return (int) mContingencyTable[NEG][POS];
  }

  /** @return the number of true positives. */
  public int falseNegatives() {
    return (int) mContingencyTable[POS][NEG];
  }

  /** @return the number of correctly classified instances. */
  public int correct() {
    return truePositives() + trueNegatives();
  }

  /** @return the number of incorrectly classified instances. */
  public int incorrect() {
    return falsePositives() + falseNegatives();
  }

  /** @return the total number of classified instances. */
  public int total() {
    return correct() + incorrect();
  }

  /** @return the fraction of incorrectly classified instances. */
  public double errorRate() {
    return (double) incorrect() / (incorrect() + correct());
  }

  /** @return the fraction of correctly classified instances. */
  public double accuracy() {
    return (double) correct() / (incorrect() + correct());
  }
}
