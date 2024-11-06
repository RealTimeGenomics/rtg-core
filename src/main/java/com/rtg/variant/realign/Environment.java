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
package com.rtg.variant.realign;


/**
 * Stores information about a read, a template and their relationship.
 *
 */
public interface Environment {

  /**
   * Get the probability of the specified nucleotide
   * at the given index in the read. These probabilities will typically
   * be computed from quality scores.
   * @param index position in read (0 based).
   * @return the raw probability (not its log).
   */
  double quality(int index);

  /**
   * Get the base at specified position in the read.
   * @param index position in read (0 based).
   * @return the nucleotide (0 = N).
   */
  byte read(int index);

  /**
   * Get the base at the specified relative position in the template.
   * @param index position in template (0 based from the expected start of the read).
   * @return the nucleotide (0 = N).
   */
  byte template(int index);

  /**
   * Convert a relative position on the template to its real position.
   * This does NOT take into account the actual alignment of the read,
   * or any CG gaps.  It simply assumes that the read is laid along the
   * template one-for-one.
   * @param index zero means the nominal 'start' of the read on the template
   *   (which will be the end of the read for a reversed read).
   * @return the real template position (may be off the template).
   */
  int absoluteTemplatePosition(int index);

  /**
   * Get a position which is guaranteed to be greater than the position of the last non-null nucleotide in the read
   * (0 based).
   * @return the length of the read.
   */
  int readLength();

  /**
   * @return the length of the template.
   */
  int templateLength();

  /**
   * Get the maximum shift permitted when doing alignment.
   * @return the maximum shift permitted when doing alignment.
   */
  int maxShift();
}
