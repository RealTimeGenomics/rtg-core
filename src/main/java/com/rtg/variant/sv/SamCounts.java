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

/**
 * Records counts of SAM records at particular coordinates.
 *
 */
public interface SamCounts {

  /**
   * Add a new record to the underlying set at position given by index.
   * @param index position used.
   */
  void increment(final int index);

  /**
   * Get the number of records at the index.
   * @param base position for window (in forward case this is at right end of
   * window).
   * @param index location being fetched (in forward case counting from left hand end of window).
   * @return count - 0 if no set at specified location
   */
  double count(int base, int index);

  /**
   * Get the number of records at all positions between <code>index</code> and <code>index2</code> (exclusive).
   * @param base position for window.
   * @param index first location being fetched (in forward case counting from left hand end of window).
   * @param index2 last location being fetched (exclusive).
   * @return count - 0 if no set at specified location
   */
  double count(int base, int index, int index2);

  /**
   * Get the sum of count * log(count) between <code>index</code> and <code>index2</code> (exclusive).
   * @param base position for window.
   * @param index first location being fetched (in forward case counting from left hand end of window).
   * @param index2 last location being fetched (exclusive).
   * @return count - 0 if no set at specified location
   */
  double sumLn(int base, int index, int index2);

  /**
   * Reset all sets at all indexes and ensure that the array
   * is able to accommodate length entries.
   * @param length to be accommodated from now on.
   * @param bufferSize size of active view
   */
  void reset(final int length, int bufferSize);

  /**
   * Sets the current minimum position that can be used
   * @param offset the current minimum position.
   */
  void flushTo(int offset);

  /**
   * The current effective length of the underlying index (used for checking index values etc.).
   * May be less than the actual backing array.
   * @return the effective length.
   */
  int length();

}
