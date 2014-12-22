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
