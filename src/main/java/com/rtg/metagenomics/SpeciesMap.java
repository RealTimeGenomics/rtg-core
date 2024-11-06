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

