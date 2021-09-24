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
package com.rtg.index.hash.ngs;

import com.rtg.index.IndexSet;

/**
 * Does the actions for each window when scanning reads.
 */
public class ReadCallImplementation implements ReadCall {
  private final IndexSet mIndexes;

  /**
   * @param indexes array of indexes selected by the hash function.
   */
  public ReadCallImplementation(final IndexSet indexes) {
    mIndexes = indexes;
  }

  @Override
  public void readCall(final int id, final long hash, final int index) {
    //Diagnostic.developerLog("add  index=" + index + " hash=" + com.rtg.util.Utils.toBitsSep(hash) + " " + hash);
    mIndexes.get(index).add(hash, id);
  }
}

