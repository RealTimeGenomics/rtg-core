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

import java.util.ArrayList;
import java.util.Arrays;

import com.rtg.metagenomics.matrix.Vector;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Information required to subset data structures belonging to blocks in the species module
 */
public class BlockMapping extends IntegralAbstract {
  /** the id of the lowest genome in the block to which the genome belongs */
  private final int[] mA;
  /** the ids of the genomes in each block */
  private final int[][] mBlocks;
  /** the position of a genome within its block */
  private final int[] mC;
  /** the block that a genome belongs to */
  private final int[] mD;

  /**
   * @param frags array of all fragments
   * @param n number of genomes
   */
  public BlockMapping(Frag[] frags, int n) {
    this(constructA(frags, n));
    assert globalIntegrity();
  }

  private BlockMapping(int[] a) {
    mA = a;
    final int n = mA.length;
    //count the number of blocks and their size
    final int[] c = new int[n];
    int bn = 0;
    final int[] cnt = new int[n];
    for (int i = 0; i < n; ++i) {
      final int ai = a[i];
      if (ai == i) {
        ++bn;
      }
      if (ai >= 0) {
        cnt[a[ai]]++;
      }
    }
    final int[][] blocks = new int[bn][];
    final int[] next = new int[bn];
    final int[] block = new int[n];
    Arrays.fill(block, -1);
    //fill in entries for each block
    int b = 0;
    for (int i = 0; i < n; ++i) {
      final int ai = a[i];
      if (cnt[i] > 0) {
        assert ai == i;
        final int numInBlock = cnt[i];
        blocks[b] = new int[numInBlock];
        block[i] = b;
        ++b;
      }
      if (ai >= 0) {
        final int bx = block[ai];
        blocks[bx][next[bx]] = i;
        c[i] = next[bx];
        block[i] = bx;
        next[bx]++;
      }

    }
    mBlocks = blocks;
    mC = c;
    mD = block;
  }

  static int[] constructA(Frag[] frags, int n) {
    final int[] a = new int[n];
    Arrays.fill(a, -1);
    for (final Frag frag : frags) {
      frag.stats(a);
    }
    for (int i = 0; i < a.length; ++i) {
      final int ai = a[i];
      if (ai >= 0 && a[ai] < ai) {
        a[i] = a[ai];
      }
    }
    return a;
  }

  private static long[] subsetGenomeLengths(int[] block, long[] genomeLengths) {
    final long[] lengths = new long[block.length];
    for (int i = 0; i < block.length; ++i) {
      lengths[i] = genomeLengths[block[i]];
    }
    return lengths;
  }

  String statistics() {
    final StringBuilder sb = new StringBuilder();
    final int n = mA.length;
    sb.append("Total genomes:").append(n).append(LS);
    sb.append("Total blocks:").append(mBlocks.length).append(LS);

    int noHit = 0;
    final int[] cnt = new int[n];
    for (final int ai : mA) {
      if (ai == -1) {
        ++noHit;
      } else {
        cnt[ai]++;
      }
    }
    final int[] cc = new int[n + 1];
    for (int i = 0; i < n; ++i) {
      cc[cnt[i]]++;
    }
    sb.append("No hits:").append(noHit).append(LS);
    for (int i = 1; i <= n; ++i) {
      final int cci = cc[i];
      if (cci > 0) {
        sb.append("block size:").append(i).append(" times:").append(cci).append(LS);
      }
    }
    return sb.toString();
  }

  /**
   * @param genomeId the global internal genome id of a species
   * @return the block to which the genome has been assigned
   */
  public int getBlockId(int genomeId) {
    return mD[genomeId];
  }

  /**
   * @param genomeId the global internal genome id of a species
   * @return the block local genome id of a species
   */
  public int getLocalId(int genomeId) {
    return mC[genomeId];
  }

  /**
   * Creates the sub blocks to perform work on.
   * @param params species parameters
   * @param frags array of all fragments
   * @param speciesMap index of all sequence species names and ids
   * @param genomeLengths length of all reference sequences
   * @return the information about each block
   */
  public BlockInfo[] subBlock(SpeciesParams params, Frag[] frags, SpeciesMap speciesMap, long[] genomeLengths) {
    final ArrayList<ArrayList<Frag>> moddedFrags = new ArrayList<>(mBlocks.length);
    for (int k = 0; k < mBlocks.length; ++k) {
      moddedFrags.add(new ArrayList<Frag>());
    }
    for (final Frag f : frags) {
      moddedFrags.get(f.subBlock(mD)).add(f.subFrag(mC));
    }
    final BlockInfo[] blockInfos = new BlockInfo[mBlocks.length];
    for (int b = 0; b < mBlocks.length; ++b) {
      final ArrayList<Frag> subFrags = moddedFrags.get(b);
      blockInfos[b] = new BlockInfo(b, speciesMap, subFrags.toArray(new Frag[subFrags.size()]), speciesMap.subset(mBlocks[b]), subsetGenomeLengths(mBlocks[b], genomeLengths), params.verbose());
    }
    return blockInfos;
  }

  /**
   * Incorporates results from <code>subResults</code> into <code>dest</code>. Translating genome indexes as required
   * @param dest destination for results
   * @param subResults results for a single block
   * @param b the block number corresponding to <code>subResults</code>
   */
  public void mergeResults(BlockResult dest, BlockResult subResults, int b) {

    // R values are just poked directly in after conversion from local to global id
    final Vector subR = subResults.getR();
    for (int i = 0; i < mBlocks[b].length; ++i) {
      dest.setR(mBlocks[b][i], subR.get(i));
    }

    // Variance log is summed, in global id space
    final Vector subStd = subResults.getVarianceLog();
    final Vector subLikelihoods = subResults.getLikelihoods();
    for (int i = 0; i < dest.getVarianceLog().size(); ++i) {
      dest.setVariance(i, dest.getVarianceLog().get(i) + subStd.get(i));
      dest.setLikelihood(i, dest.getLikelihoods().get(i) + subLikelihoods.get(i));
    }
  }

  public int[] getA() {
    return mA;
  }

  public int[][] getBlocks() {
    return mBlocks;
  }

  public int[] getC() {
    return mC;
  }

  @Override
  public boolean globalIntegrity() {
    for (int i = 0; i < mA.length; ++i) {
      Exam.assertTrue(mA[i] <= i);
      Exam.assertTrue(mA[i] == -1 || mA[mA[i]] == mA[i]);
    }

    Exam.assertEquals(mA.length, mC.length);
    for (final int[] block : mBlocks) {
      for (int i = 0; i < block.length; ++i) {
        final int x = block[i];
        Exam.assertEquals(block[0], mA[x]);
        Exam.assertEquals(i, mC[x]);
      }
    }
    return integrity();
  }

  @Override
  public boolean integrity() {
    return true;
  }
}
