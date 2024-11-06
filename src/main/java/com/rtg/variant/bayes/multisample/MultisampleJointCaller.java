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

package com.rtg.variant.bayes.multisample;


import java.io.Closeable;
import java.util.List;

import com.rtg.variant.Variant;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.snp.HypothesesPrior;

/**
 * Primary interface for multi-sample callers.
 */
public interface MultisampleJointCaller extends Closeable {

  /**
   * Make the actual call.
   * @param <D> the type of description.
   * @param <T> the type of the haploid-diploid hypotheses.
   *
   * @param templateName name of template sequence.
   * @param position zero-based position of call.
   * @param endPosition zero-based end position of call.
   * @param ref the reference genome, 0=N, 1=A, 2=C, 3=G, 4=T.
   * @param models individual models used to construct joint model.
   * @param hypotheses hypotheses containing priors
   * @return Variant object, or null if call is not to be retained
   */
  <D extends Description, T extends HypothesesPrior<D>> Variant makeCall(String templateName, int position, int endPosition, byte[] ref, List<ModelInterface<?>> models, HaploidDiploidHypotheses<T> hypotheses);

  /**
   * Called at the end of each input sequence.
   */
  void endOfSequence();

}
