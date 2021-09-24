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

import com.rtg.util.ContingencyTable;

/**
 * Evaluate a classifier with respect to a test dataset.
 */
public class SimpleEvaluation extends ContingencyTable {

  /**
   * Evaluate the classifier on the supplied dataset
   * @param classifier the classifier to evaluate
   * @param dataset the dataset to evaluate on
   */
  public void evaluate(PredictClassifier classifier, Dataset dataset) {
    for (Instance instance : dataset.getInstances()) {
      final int actual = instance.isPositive() ? POS : NEG;
      final int predicted = classifier.predict(instance.instance()) > 0.5 ? POS : NEG;
      add(actual, predicted, instance.weight());
    }
  }

}
