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
