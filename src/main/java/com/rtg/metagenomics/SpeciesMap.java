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
package com.rtg.metagenomics;

import java.util.HashMap;
import java.util.Map;

/**
 * Mapping from species sequence names to incremental integer identifier for each species.
 */
public final class SpeciesMap extends HashMap<Integer, Integer> {

  private int mNextId = 0;
  private final HashMap<Integer, Integer> mBack = new HashMap<>();

  /**
   * @param taxonId name of the species.
   * @return the unique identifier for the species.
   */
  int id(final Integer taxonId) {
    final Integer i = this.get(taxonId);
    if (i == null) {
      this.put(taxonId, mNextId);
      mBack.put(mNextId, taxonId);
      return mNextId++;
    } else {
      return i;
    }
  }

  /**
   * Gets the taxon id corresponding to the given unique identifier
   * @param id the unique identifier
   * @return the taxon id
   */
  int taxonId(int id) {
    return mBack.get(id);
  }

  int[] taxonIds() {
    final int[] n = new int[mBack.size()];
    for (final Map.Entry<Integer, Integer> e : mBack.entrySet()) {
      n[e.getKey()] = e.getValue();
    }
    return n;
  }

  SpeciesMap subset(int... genomes) {
    final SpeciesMap ret = new SpeciesMap();
    for (int g : genomes) {
      ret.id(mBack.get(g));
    }
    return ret;
  }
}

