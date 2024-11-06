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
