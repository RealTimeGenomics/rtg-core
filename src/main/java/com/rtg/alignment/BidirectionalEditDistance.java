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
package com.rtg.alignment;

/**
 * Interface for calculating an edit distance measure between two strings
 * in either forward or reverse directions
 */
public interface BidirectionalEditDistance {

  /**
   * Calculate a global alignment between an array and a specified portion
   * of another array. Higher weights, cause greater penalty for a particular
   * type of operation. There is an initial penalty for starting an indel
   * run, thereafter they suffer only ordinary penalty.
   * Ns added on either side of the template suffer an ordinary penalty.
   * @param read the read sequence
   * @param rlen length of read (assumed that given read is at least this long)
   * @param template the template sequence
   * @param zeroBasedStart the start position
   * @param rc true if a reverse complement alignment should be done
   * @param maxScore maximum score for an alignment to allow early termination
   * @param maxShift maximum allowed start position alignment shift
   * @param cgLeft true if this is the left arm (used only in CG case so far)
   * @return packed alignment information
   */
  int[] calculateEditDistance(byte[] read, int rlen, byte[] template, int zeroBasedStart, boolean rc, int maxScore, int maxShift, boolean cgLeft);

  /**
   * Log statistics from the alignment run.
   */
  void logStats();

}
