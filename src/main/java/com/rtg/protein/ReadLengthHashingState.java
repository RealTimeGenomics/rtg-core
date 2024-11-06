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
package com.rtg.protein;

import com.rtg.index.IndexSet;
import com.rtg.index.hash.ngs.NgsHashFunction;
import com.rtg.index.hash.ngs.ProteinNgsHashLoop;
import com.rtg.index.hash.ngs.TemplateCall;

/**
 * Various state associated with each read length for protein mapping.
 */
public class ReadLengthHashingState {
  private final IndexSet mIndexes;
  private final ProteinNgsHashLoop mHashLoop;
  private final NgsHashFunction mHashFunction;
  private final TemplateCall mTemplateCallImplementation;

  /**
   * @param indexSet indexes for this read length
   * @param loop hash loop for the read length
   * @param hashFunction hash function for this read length
   * @param tci template call for this read length
   */
  public ReadLengthHashingState(IndexSet indexSet, ProteinNgsHashLoop loop, NgsHashFunction hashFunction, TemplateCall tci) {
    mIndexes = indexSet;
    mHashLoop = loop;
    mHashFunction = hashFunction;
    mTemplateCallImplementation = tci;
  }

  /**
   * @return the hash function
   */
  public NgsHashFunction getHashFunction() {
    return mHashFunction;
  }

  /**
   * @return the hash loop
   */
  public ProteinNgsHashLoop getHashLoop() {
    return mHashLoop;
  }

  /**
   * @return the indexes for the read
   */
  public IndexSet getIndexes() {
    return mIndexes;
  }

  /**
   * @return the template call implementation for the chunk
   */
  public TemplateCall getTemplateCallImplementation() {
    return mTemplateCallImplementation;
  }

}
