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
