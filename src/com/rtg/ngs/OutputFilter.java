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
    public OutputProcessor makeProcessor(final NgsParams params, MapStatistics stats) throws IOException {
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

