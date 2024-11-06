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

package com.rtg.variant;

import java.util.Collection;

import com.rtg.relation.GenomeRelationships;

/**
 * Enumeration with possible genome connectivity assignments for a pedigree.
 */
public enum GenomeConnectivity {
  /** Genomes are poorly connected */
  SPARSE,
  /** Genomes are well connected */
  DENSE;

  /**
   * Determines the connectivity of the given genome relationships based on the genomes that are available to analyse.
   * @param genomes list of genome samples for which there is data to analyse.
   * @param pedigree relationships between the genomes.
   * @return the connectivity of the genomes.
   */
  public static GenomeConnectivity getConnectivity(Collection<String> genomes, GenomeRelationships pedigree) {
    // TODO: workout how connected the pedigree is
    // if 2 or more groups/individuals then assume sparse ???
    return pedigree.numberOfDisconnectedGroups(genomes) <= 2 ? DENSE : SPARSE;
    //return SPARSE;
  }
}
