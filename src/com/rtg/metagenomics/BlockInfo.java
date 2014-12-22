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

import java.util.Map.Entry;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.metagenomics.matrix.Vector;
import com.rtg.metagenomics.matrix.VectorSimple;
import com.rtg.util.SortedMultiSet;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Holds information needed for a solving a species block.
 */
@TestClass("com.rtg.metagenomics.SpeciesTest")
public class BlockInfo extends IntegralAbstract {

  /** If set, output even more information than usual */
  public static final boolean VERY_VERBOSE = Boolean.getBoolean("rtg.species.very-verbose");

  private final int mBlockId;

  private final boolean mVerbose;

  private final SpeciesMap mGlobalSpeciesMap; // null if the ids here are the same as the global ids

  private final Frag[] mFrags;

  private final SpeciesMap mSpeciesMap;

  private final long[] mGenomeLengths;

  private final Vector mGenomeLengthsVector;

  private final int[] mTaxonIds;

  private final int mN;

  private final int mTotalReads;

  /**
   * @param blockId unique id for this block.
   * @param globalSpeciesMap translates global species id to taxon id
   * @param frags fragments for this block
   * @param speciesMap translates block-specific species id to taxon id
   * @param genomeLengths lengths of species for this block, indexed by species id
   * @param verbose if true, output extra information
   */
  BlockInfo(int blockId, SpeciesMap globalSpeciesMap, Frag[] frags, SpeciesMap speciesMap, final long[] genomeLengths, boolean verbose) {
    mBlockId = blockId;
    mVerbose = verbose;
    mGlobalSpeciesMap = globalSpeciesMap;
    mFrags = frags;
    mSpeciesMap = speciesMap;
    mGenomeLengths = genomeLengths;
    mN = genomeLengths.length;
    mGenomeLengthsVector = new VectorSimple(mN);
    for (int i = 0; i < genomeLengths.length; i++) {
      mGenomeLengthsVector.incr(i, genomeLengths[i]);
    }
    mTaxonIds = mSpeciesMap == null ? null : mSpeciesMap.taxonIds();
    int reads = 0;
    for (final Frag frag : frags) {
      reads += frag.multiplicity();
    }
    mTotalReads = reads;
  }


  /**
   * Get the length of the <code>i'th</code> genome in this block.
   * @param i index.
   * @return length of the <code>i'th</code> genome in this block.
   */
  public long getGenomeLength(final int i) {
    return mGenomeLengths[i];
  }

  /**
   * Get the taxon id of the <code>i'th</code> genome in this block.
   * @param i index.
   * @return taxon id of the <code>i'th</code> genome in this block.
   */
  public int getTaxonId(int i) {
    return mTaxonIds[i];
  }

  public SpeciesMap getGlobalSpeciesMap() {
    return mGlobalSpeciesMap;
  }

  public int getN() {
    return mN;
  }

  public Vector getGenomeLengthsVector() {
    return mGenomeLengthsVector;
  }

  public boolean isVerbose() {
    return mVerbose;
  }

  public Frag[] getFrags() {
    return mFrags;
  }

  public SpeciesMap getSpeciesMap() {
    return mSpeciesMap;
  }

  @Override
  public boolean integrity() {
    return false;
  }

  /**
   * @return a human readable summary of the information in the fragments.
   */
  String fragInfo() {
    final StringBuilder sb = new StringBuilder();
    @SuppressWarnings("unchecked")
    final SortedMultiSet<Integer>[] counts = (SortedMultiSet<Integer>[]) new SortedMultiSet<?>[mTaxonIds.length];
    for (int i = 0; i < counts.length; i++) {
      counts[i] = new SortedMultiSet<>();
    }
    for (Frag mFrag : mFrags) {
      mFrag.count(counts);
    }
    sb.append("block ").append(mBlockId).append(LS);
    for (int i = 0; i < counts.length; i++) {
      sb.append(i).append(": ").append(mTaxonIds[i]).append(" {");
      final SortedMultiSet<Integer> cntSet = counts[i];
      for (final Entry<Integer, Integer> entry : cntSet.countMap().entrySet()) {
        sb.append(entry.getKey()).append("#").append(entry.getValue()).append(", ");
      }
      sb.append("}").append(LS);
    }
    return sb.toString();
  }

  @Override
  public String toString() {
    return "BlockInfo{"
        + "mBlockId=" + mBlockId
        + ", mVerbose=" + mVerbose
        + ", mN=" + mN
        + ", mTotalReads=" + mTotalReads
        + '}';
  }


  /**
   * @return unique identifier for this block.
   */
  public int id() {
    return mBlockId;
  }
}
