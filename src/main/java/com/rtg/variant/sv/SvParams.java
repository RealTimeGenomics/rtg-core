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

package com.rtg.variant.sv;

import java.util.Collections;
import java.util.Map;

import com.rtg.sam.SingleMappedParams;

/**
 * Common Structural Variant Parameters
 */
public abstract class SvParams extends SingleMappedParams {

  /**
   * A builder class for <code>SVParams</code>.
   * @param <B> builder type
   */
  public abstract static class SvParamsBuilder<B extends SvParamsBuilder<B>> extends SingleMappedParamsBuilder<B> {
    protected Map<String, String> mReadGroupLabels;
    protected Map<String, ReadGroupStats> mReadGroupStatistics;

    /**
     * Set the read group re-labeling map
     * @param labels the relabeling map
     * @return this builder, so calls can be chained.
     */
    public B readGroupLabels(Map<String, String> labels) {
      mReadGroupLabels = labels;
      return self();
    }

    /**
     * Set the read group stats
     * @param stats a pre build read group stats map
     * @return this builder, so calls can be chained.
     */
    public B readGroupStatistics(Map<String, ReadGroupStats> stats) {
      mReadGroupStatistics = stats;
      return self();
    }
  }

  private final Map<String, String> mReadGroupLabels;
  private final Map<String, ReadGroupStats> mReadGroupStatistics;

  /**
   * @param builder the builder object.
   */
  protected SvParams(SvParamsBuilder<?> builder) {
    super(builder);
    mReadGroupLabels = builder.mReadGroupLabels == null ? null : Collections.unmodifiableMap(builder.mReadGroupLabels);
    mReadGroupStatistics = builder.mReadGroupStatistics == null ? null : Collections.unmodifiableMap(builder.mReadGroupStatistics);
  }

  /**
   * Get per read-group statistics.
   * @return mapping to read group statistics.
   */
  public Map<String, ReadGroupStats> readGroupStatistics() {
    return mReadGroupStatistics;
  }

  /**
   * Get the read-group labelling map.
   * @return mapping from input to output read group names.
   */
  public Map<String, String> readGroupLabels() {
    return mReadGroupLabels;
  }

}
