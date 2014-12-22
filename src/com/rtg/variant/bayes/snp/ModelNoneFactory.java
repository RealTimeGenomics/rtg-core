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

package com.rtg.variant.bayes.snp;

import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.ModelFactory;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.ModelNone;

/**
 */
public class ModelNoneFactory implements ModelFactory<Description, HypothesesNone<Description>> {

  @Override
  public ModelInterface<Description> make(int ref) {
    return ModelNone.SINGLETON;
  }

  @Override
  public HypothesesNone<Description> defaultHypotheses(int ref) {
    return HypothesesNone.SINGLETON;
  }

}
