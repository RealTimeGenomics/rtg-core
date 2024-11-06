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
package com.rtg.position;

import com.rtg.index.Add;
import com.rtg.index.IndexUtils;
import com.rtg.index.hash.ExactHashFunction;
import com.rtg.index.hash.HashFunction;
import com.rtg.index.hash.HashLoop;
import com.rtg.index.hash.InExactHashFunction;
import com.rtg.index.params.ParamsUtils;
import com.rtg.launcher.BuildParams;
import com.rtg.position.output.PositionParams;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.SlimException;

/**
 * Takes subjects, queries and writes to "out" in output directory.
 * Also accepts the various other parameters required.
 * Used to be Position Task
 */
public final class PositionUtils {

  private PositionUtils() { }

  /**
   * Get  a string with a  human readable description of the memory usage.
   * It conforms to the usage that each line starts with the word "Memory" followed by
   * a unique identifying string that doesn't contain spaces followed by the
   * number of bytes. All entries are to be disjoint (so no totals are to be included).
   * @param buildSearchParams parameters to be used to create <code>BuildSearch</code>
   * @return a string with a  human readable description of the memory usage.
   */
  public static String memToString(final PositionParams buildSearchParams) {
    final StringBuilder sb = new StringBuilder();
    memToString(sb, buildSearchParams);
    return sb.toString();
  }

  /**
   * Get  a string with a  human readable description of the memory usage.
   * It conforms to the usage that each line starts with the word "Memory" followed by
   * a unique identifying string that doesn't contain spaces followed by the
   * number of bytes. All entries are to be disjoint (so no totals are to be included).
   * @param buildSearchParams parameters to be used to create <code>BuildSearch</code>
   * @param sb where to put the result.
   */
  private static void memToString(final StringBuilder sb, final PositionParams buildSearchParams) {
    sb.append(ParamsUtils.memToString("Shared_buffer", buildSearchParams.bufferLength()));
    sb.append(IndexUtils.memToString(buildSearchParams.indexParams()));
  }


  /**
   * Construct a <code>HashLoop</code> to use when doing paired end processing.
   * @param index the index of hashes.
   * @param buildFirst the parameters for one arm of the pair.
   * @param buildSecond  the parameters for the other arm of the pair.
   * @param winBits number of bits required to encode a window
   * @return the <code>HashLoop</code>.
   */
  public static HashLoop makePairedBuild(final Add index, final BuildParams buildFirst, final BuildParams buildSecond, int winBits) {
    final int mxs = Math.max(maxMatches(buildFirst), maxMatches(buildSecond));
    return makeBuild(index, buildFirst, winBits, mxs, true);
  }

  /**
   * Construct a <code>HashLoop</code> to use when doing single end processing.
   * @param index the index of hashes.
   * @param buildParams the parameters.
   * @param winBits number of bits required to encode a window
   * @return the <code>HashLoop</code>.
   */
  public static HashLoop makeBuild(final Add index, final BuildParams buildParams, int winBits) {
    final int mxs = maxMatches(buildParams);
    return makeBuild(index, buildParams, winBits, mxs, false);
  }

  private static HashLoop makeBuild(final Add index, final BuildParams buildParams, final int winBits, final int mxs, final boolean pairMode) {
    final HashLoop subjectHashLoop;
    if (winBits > 64) {
      final HashFunction function = new InExactHashFunction(buildParams.windowSize());
      subjectHashLoop = new BuildResetHashLoop(buildParams.windowSize(), buildParams.stepSize(), function, index, mxs, pairMode);
    } else {
      final ExactHashFunction exf = new ExactHashFunction(buildParams);
      subjectHashLoop = new BuildIncrementalHashLoop(buildParams.stepSize(), exf, index, mxs, pairMode);
    }
    return subjectHashLoop;
  }

  /**
   * Compute the number of different match positions that can occur in a build sequence.
   * @param buildParams parameters.
   * @return the number of possible match positions (&ge; 0).
   * @throws IllegalArgumentException if result is too large.
   */
  public static int maxMatches(final BuildParams buildParams) {
    final int stepSize = buildParams.stepSize();
    final int windowSize = buildParams.windowSize();
    final long maxLength = buildParams.sequences().maxLength();
    return maxMatches(windowSize, stepSize, maxLength);
  }

  /**
   * Compute the number of different match positions that can occur in a build sequence.
   * @param windowSize window size parameter
   * @param stepSize step size parameter
   * @param maxLength length of longest sequence
   * @return the number of possible match positions (&ge; 0).
   * @throws IllegalArgumentException if result is too large.
   */
  public static int maxMatches(int windowSize, int stepSize, long maxLength) {
    final long mxs;
    if (maxLength < windowSize) {
      mxs = 0;
    } else {
      mxs = (maxLength - windowSize) / stepSize + 1;
    }
    assert mxs >= 0 : mxs;
    if (mxs > Integer.MAX_VALUE) {
      throw new IllegalArgumentException("Unable to store identifier bits. Input data too long. max sequence length=" + maxLength + " step=" + stepSize + " window=" + windowSize);
    }
    return (int) mxs;
  }

  /**
   * Make a buffer long enough to be used for both build and search.
   * @param maxSequence the length of the longest sequence.
   * @return an array long enough to be used as a buffer.
   * @throws SlimException if <code>maxSequence</code> is too long.
   */
  public static byte[] makeBuffer(final long maxSequence) {
    if (maxSequence > Integer.MAX_VALUE) {
      Diagnostic.error(ErrorType.SEQUENCE_TOO_LONG, String.valueOf(maxSequence));
      throw new SlimException();
    }
    return new byte[(int) maxSequence];
  }
}
