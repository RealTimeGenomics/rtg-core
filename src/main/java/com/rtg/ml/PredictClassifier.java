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

import java.io.DataOutputStream;
import java.io.IOException;

/**
 * Two-class instance classification.
 *
 */
public interface PredictClassifier {

  /**
   * Save the classifier to given stream
   * @param dos the stream to save to
   * @param data attribute set to help encode and decode values
   * @throws IOException if an IO error occurs
   */
  void save(DataOutputStream dos, Dataset data) throws IOException;

  /**
   * Return the classifier probability that the instance is an exemplar of the positive class.
   * @param instance the instance
   * @return probability of positive
   */
  double predict(double[] instance);

  /**
   * Get a human readable representation of the classifier
   * @param out where to send the output
   * @param indent an indentation level that should be used for each output line
   * @param data set of attributes for encoding/decoding values
   * @return the output StringBuilder, for convenience.
   */
  StringBuilder toString(StringBuilder out, String indent, Dataset data);

}
