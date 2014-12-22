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
