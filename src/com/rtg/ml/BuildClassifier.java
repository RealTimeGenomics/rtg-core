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

import java.util.Properties;

/**
 * Build a two-class classifier.
 *
 */
public interface BuildClassifier {

  /**
   * Train the classifier on the specified dataset.
   * @param dataset training data
   */
  void build(Dataset dataset);

  /**
   * Return the classifier that has been built.
   * @return prediction classifier.
   */
  PredictClassifier getClassifier();

  /**
   * Sets additional properties (if any) associated with this classifier.
   * @param properties any properties that alter how the model gets built.
   */
  void setProperties(Properties properties);

}
