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
package com.rtg.index;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.rtg.index.params.CreateParams;
import com.rtg.util.diagnostic.Diagnostic;

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
    Diagnostic.userLog("Creating blacklist index");
    final IndexCompressed blacklistIndex = new IndexCompressed(createParams, new UnfilteredFilterMethod(), numberThreads);
    for (long hash : blacklist) {
      blacklistIndex.add(hash, 0);
    }
    blacklistIndex.freeze();
    for (long hash : blacklist) {
      blacklistIndex.add(hash, 0);
    }
    blacklistIndex.freeze();
    Diagnostic.userLog("Blacklist index constructed");
    mBlacklist = blacklistIndex;
  }

  private BlacklistFilterMethod(Index blacklist) {
    mBlacklist = blacklist;
  }

  @Override
  public IndexFilterMethod threadClone() {
    return new BlacklistFilterMethod(mBlacklist);
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

  @Override
  public String toString() {
    return "Blacklist";
  }
}
