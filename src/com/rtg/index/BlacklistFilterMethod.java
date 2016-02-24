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
package com.rtg.index;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.rtg.index.params.CreateParams;

/**
 * Filter which uses a blacklist
 */
public class BlacklistFilterMethod implements IndexFilterMethod {

  private final Index mBlacklist;

  /**
   * @param blacklist the list of blacklisted hashes
   * @param hashBits how many bits of each long is used to encode hash
   * @param numberThreads number of threads for index sorting
   */
  public BlacklistFilterMethod(List<Long> blacklist, int hashBits, int numberThreads) {
    final CreateParams createParams = new CreateParams.CreateParamsBuilder().valueBits(0).compressHashes(true).createBitVector(true).windowBits(hashBits).size(blacklist.size()).hashBits(hashBits).create();
    final IndexCompressed blacklistIndex = new IndexCompressed(createParams, new UnfilteredFilterMethod(), numberThreads);
    for (long hash : blacklist) {
      blacklistIndex.add(hash, 0);
    }
    blacklistIndex.freeze();
    for (long hash : blacklist) {
      blacklistIndex.add(hash, 0);
    }
    blacklistIndex.freeze();
    mBlacklist = blacklistIndex;
  }

  /**
   *
   * @param sdfDir directory of reference SDF
   * @param wordSize kmer size for hashes
   * @param threshold kmer occurence count at which they get placed in blacklist
   * @param numberThreads number of threads for index sorting
   * @return the filter
   * @throws IOException if an IO error occurs
   */
  public static BlacklistFilterMethod loadBlacklist(File sdfDir, int wordSize, int threshold, int numberThreads) throws IOException {
    final List<Long> blacklist = HashBlacklist.loadBlacklist(sdfDir, wordSize, threshold);
    return new BlacklistFilterMethod(blacklist, HashBlacklist.hashBits(wordSize), numberThreads);
  }

  @Override
  public void initialize(Index index) {
  }

  @Override
  public boolean keepHash(long hash, long numHits) {
    return !mBlacklist.contains(hash);
  }
}
