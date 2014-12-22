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

package com.rtg.variant.eval;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.vcf.VcfRecord;

/**
 * Class template for value extractor.
 */
@TestClass("com.rtg.variant.eval.RocScoreFieldTest")
public abstract class RocSortValueExtractor {

  abstract double getSortValue(VcfRecord rec, int sampleNo);

  /** Dummy null extractor for testing purposes */
  public static final RocSortValueExtractor NULL_EXTRACTOR = new RocSortValueExtractor() {
    @Override
    double getSortValue(VcfRecord rec, int sampleNo) {
      return 0;
    }
  };
}
