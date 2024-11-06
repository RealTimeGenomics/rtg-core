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
package com.rtg.index.hash.ngs;

/**
 * Creates new copies of <code>SplitHashLoop</code> when requested.
 * Also supplies the number of windows used (needed so that the indexes can be created before
 * the actual <code>SplitHashLoop</code> is available.
 */
public interface HashFunctionFactory {

  /**
   * Get the number of windows that will be used by the <code>SplitHashLoop</code>.
   * @return the number of windows that will be used by the <code>SplitHashLoop</code>.
   */
  int numberWindows();

  /**
   * Get the number of bits recorded in each hash window.
   * @return the number of bits recorded in each hash window.
   */
  int hashBits();

  /**
   * Get the number of bits needed to uniquely represent a window (the hash may lose information).
   * @return the number of bits needed to uniquely represent a window.
   */
  int windowBits();

  /**
   * Get the window size
   * @return the size
   */
  int windowSize();

  /**
   * Create a new instance of a <code>HashFunction</code>.
   * @param readCall will do processing of read sequences (during build).
   * @param templateCall will do processing of template sequences (during search).
   * @return a new instance of a <code>HashFunction</code>.
   */
  NgsHashFunction create(ReadCall readCall, TemplateCall templateCall);
}

