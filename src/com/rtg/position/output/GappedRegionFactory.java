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
package com.rtg.position.output;

import com.rtg.util.array.ImmutableIntArray;

/**
 * factory for gapped region objects
 * @param <G> Type of region object to being created
 */
public interface GappedRegionFactory<G extends AbstractGappedRegion<G>> {

  /**
   * Create a gapped region appropriate for this output type.
   * If regions are not appropriate then the default is to throw an <code>IllegalOperationException</code>.
   * @param id an identifier unique within a particular <code>GappedOutput</code> object.
   * @param distribution probability distribution for gap probabilities.
   * @param params parameters for current job.
   * @param buildLengths build lengths.
   * @return an appropriate gapped region (never null).
   */
  G region(final int id, final GapScorer distribution, final PositionParams params, final ImmutableIntArray buildLengths);
}
