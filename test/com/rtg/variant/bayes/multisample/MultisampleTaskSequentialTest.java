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

import com.rtg.variant.ThreadingEnvironment;
import com.rtg.variant.VariantParamsBuilder;

/**
 */
public class MultisampleTaskSequentialTest extends MultisampleTaskTest {

  @Override
  protected VariantParamsBuilder getBuilder() {
    final VariantParamsBuilder builder = super.getBuilder();
    builder.threadingEnvironment(ThreadingEnvironment.SINGLE);
    return builder;
  }

  @Override
  protected String[] getArguments() {
    return new String[] {"--Xthreading-env", "single"};
  }

}
