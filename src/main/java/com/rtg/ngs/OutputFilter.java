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
package com.rtg.ngs;


import java.io.IOException;
import java.io.Serializable;
import java.util.Locale;

import com.rtg.index.hash.ngs.OutputProcessor;
import com.rtg.protein.AllHitsProteinOutputProcessor;
import com.rtg.protein.TopEqualProteinOutputProcessor;
import com.rtg.protein.TopNProteinOutputProcessor;
import com.rtg.util.EnumHelper;
import com.rtg.util.PseudoEnum;

/**
 */
public abstract class OutputFilter implements PseudoEnum, Serializable {

  private static int sCounter = -1;

  /** No filter all the results are written. */
  public static final OutputFilter NONE = new OutputFilter(++sCounter, "NONE") {
    @Override
    public OutputProcessor makeProcessor(final NgsParams params, MapStatistics stats) throws IOException {
      return new DefaultOutputProcessorSynch(params);
    }
  };


  /** Produces SAM output for paired end data, only mated output. */
  public static final OutputFilter PAIRED_END = new OutputFilter(++sCounter, "PAIRED_END") {
    @Override
    public OutputProcessor makeProcessor(final NgsParams params, MapStatistics stats) throws IOException {
      return new TopNPairedEndOutputProcessorSync(params, stats, false, false);
    }
  };

  /** Produces SAM output for paired end data, both mated and unmated output. */
  public static final OutputFilter TOPN_PAIRED_END = new OutputFilter(++sCounter, "TOPN_PAIRED_END") {
    @Override
    public OutputProcessor makeProcessor(final NgsParams params, MapStatistics stats) throws IOException {
      return new TopNPairedEndOutputProcessorSync(params, stats, params.outputParams().outputUnmated(), params.outputParams().outputUnmapped());
    }
  };

  /** No filter all the results are written. */
  public static final OutputFilter PROTEIN_ALL_HITS = new OutputFilter(++sCounter, "PROTEIN_ALL_HITS") {
    @Override
    public OutputProcessor makeProcessor(final NgsParams params, final MapStatistics stats) throws IOException {
      return new AllHitsProteinOutputProcessor(params, stats);
    }
  };

  /** Output the N results equal to the top score. */
  public static final OutputFilter PROTEIN_TOPEQUAL = new OutputFilter(++sCounter, "PROTEIN_TOPEQUAL") {
    @Override
    public OutputProcessor makeProcessor(final NgsParams params, MapStatistics stats) throws IOException {
      return new TopEqualProteinOutputProcessor(params, stats);
    }
  };

  /** Output the best N results. */
  public static final OutputFilter PROTEIN_TOPN = new OutputFilter(++sCounter, "PROTEIN_TOPN") {
    @Override
    public OutputProcessor makeProcessor(final NgsParams params, MapStatistics stats) throws IOException {
      if (params.outputParams().filter().topN() == 1) {
        return new TopEqualProteinOutputProcessor(params, stats);
      } else {
        return new TopNProteinOutputProcessor(params, stats);
      }
    }
  };

  /** Produce SAM output for single end data. */
  public static final OutputFilter SAM_SINGLE_END = new OutputFilter(++sCounter, "SAM_SINGLE_END") {
    @Override
    public OutputProcessor makeProcessor(final NgsParams params, MapStatistics stats) throws IOException {
      return new SamSingleEndOutputProcessor(params, stats, params.outputParams().outputUnmapped());
    }
  };

  /** Produce SAM unfiltered output. */
  public static final OutputFilter SAM_UNFILTERED = new OutputFilter(++sCounter, "SAM_UNFILTERED") {
    @Override
    public OutputProcessor makeProcessor(final NgsParams params, MapStatistics stats) throws IOException {
      if (params.paired()) {
        return new UnfilteredPairedEndOutputProcessor(params, stats, params.outputParams().outputUnmapped());
      }
      return new UnfilteredSingleEndOutputProcessor(params, stats, params.outputParams().outputUnmapped());
    }
  };

  /** Produce SAM no output, for testing speed of indexing / searching. */
  public static final OutputFilter NULL = new OutputFilter(++sCounter, "NULL") {
    @Override
    public OutputProcessor makeProcessor(final NgsParams params, MapStatistics stats) {
      return new NullOutputProcessor();
    }
  };

  private final int mOrdinal;
  private final String mName;

  private OutputFilter(final int ordinal, final String name) {
    mOrdinal = ordinal;
    mName = name;
  }

  @Override
  public int ordinal() {
    return mOrdinal;
  }

  @Override
  public String name() {
    return mName;
  }

  @Override
  public String toString() {
    return mName;
  }

  private static final EnumHelper<OutputFilter> HELPER = new EnumHelper<>(OutputFilter.class, new OutputFilter[] {NONE, PAIRED_END, TOPN_PAIRED_END, PROTEIN_ALL_HITS, PROTEIN_TOPEQUAL, PROTEIN_TOPN, SAM_SINGLE_END, SAM_UNFILTERED, NULL});

  /**
   * see {@link java.lang.Enum#valueOf(Class, String)}
   * @param str name of value
   * @return the enum value
   */
  public static OutputFilter valueOf(final String str) {
    return HELPER.valueOf(str.toUpperCase(Locale.getDefault()));
  }

  /**
   * @return list of enum values
   */
  public static OutputFilter[] values() {
    return HELPER.values();
  }

  @Override
  public boolean equals(final Object arg0) {
    return this == arg0;
  }

  @Override
  public int hashCode() {
    return ordinal() + 1;
  }

  /**
   * Make an output processor.
   * @param params the parameters.
   * @param stats map to put statistics into
   * @return output processor.
   * @throws IOException if the <code>OutputProcessor</code> has problems setting up.
   */
  public abstract OutputProcessor makeProcessor(NgsParams params, MapStatistics stats) throws IOException;

}

