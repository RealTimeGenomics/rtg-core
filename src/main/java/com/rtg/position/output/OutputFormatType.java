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
package com.rtg.position.output;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import com.rtg.util.array.ImmutableIntArray;

/**
 * Output format type.
 */
public abstract class OutputFormatType implements Serializable {

  private static class Segment extends OutputFormatType {
    Segment() {
      super(0, "SEGMENT");
    }
    @Override
    public
    PositionOutput output(final PositionParams params, final GapBucketsInfo bucketInfo, final ImmutableIntArray readLengths, final Appendable out, final Appendable unmappedOut, PositionWriter writer, int repeatFreq) {
      return new SegmentOutput(params, out);
    }

  }

  /** Segmented output. Adjacent or overlapping word hits are merged. */
  public static final OutputFormatType SEGMENT = new Segment();


  private static class NgsOutputFormatType extends OutputFormatType {

    NgsOutputFormatType(final int ordinal, final String name) {
      //TODO find out if should write scoring
      super(ordinal, name);
    }

    @Override
    public
    PositionOutput output(final PositionParams params, final GapBucketsInfo bucketInfo, final ImmutableIntArray readLengths, final Appendable out, final Appendable unmappedOut, PositionWriter writer, int repeatFreq) {
      final GapScorer gs = new MismatchScores(params.output().maxGap(), params.output().distribution().maxIndel(), params.build().windowSize(), params.build().stepSize());
      return new GappedOutput<>(params, GappedScoreLongRead.FACTORY, gs, readLengths, out, writer, bucketInfo);
    }

  }

  /** Long read output format type for calling position from Ngs */
  public static final OutputFormatType NGS = new NgsOutputFormatType(1, "NGS");

  private final int mOrdinal;
  private final String mName;

  /**
   * @param name of the format
   * @param ordinal an integer id
   */
  private OutputFormatType(final int ordinal, final String name) {
    mOrdinal = ordinal;
    mName = name;
  }

  /**
   * @param params parameters.
   * @param bucketInfo bucket info
   * @param readLengths lengths of reads (only used by some of the scoring options).
   * @param out where to write the output.
   * @param unmappedOut write unmapped sequences here.
   * @param writer writer to output to.
   * @param maxHashCount max hash count.
   * @return the object that writes the output.
   */
  public abstract PositionOutput output(PositionParams params, GapBucketsInfo bucketInfo, ImmutableIntArray readLengths, Appendable out, Appendable unmappedOut, PositionWriter writer, int maxHashCount);

  private static final List<OutputFormatType> VALUES = new ArrayList<>();
  static {
    VALUES.add(SEGMENT);
    VALUES.add(NGS);
  }

  /**
   * @return list of output format types
   */
  public static List<OutputFormatType> values() {
    return new ArrayList<>(VALUES);
  }

  /**
   * @return the ordinal
   */
  public int ordinal() {
    return mOrdinal;
  }
  @Override
  public String toString() {
    return mName;
  }

  /**
   * @return name of the Output format
   */
  public String name() {
    return mName;
  }
}

