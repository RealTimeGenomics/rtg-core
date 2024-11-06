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
package com.rtg.blacklist;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Map;
import java.util.TreeMap;

import com.rtg.index.HashBlacklist;
import com.rtg.index.hash.ExactHashFunction;
import com.rtg.index.hash.HashLoop;
import com.rtg.index.hash.IncrementalHashLoop;
import com.rtg.index.params.CreateParams;
import com.rtg.launcher.BuildParams;
import com.rtg.launcher.HashingRegion;
import com.rtg.launcher.NoStatistics;
import com.rtg.launcher.ParamsTask;
import com.rtg.mode.DNA;
import com.rtg.usage.UsageMetric;
import com.rtg.util.SimpleThreadPool;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.diagnostic.OneShotTimer;

/**
 * Task for producing k-mer blacklists and a count histogram
 */
public class HashDist extends ParamsTask<HashDistParams, NoStatistics> {

  private static final String BLACKLIST_FILENAME = "blacklist";

  /**
   * @param params parameters
   * @param reportStream unused
   * @param stats unused
   * @param usageMetric unused
   */
  public HashDist(HashDistParams params, OutputStream reportStream, NoStatistics stats, UsageMetric usageMetric) {
    super(params, reportStream, stats, usageMetric);
  }

  @Override
  protected void exec() throws IOException {
    hashCount(mParams);
  }

  static void hashCount(final HashDistParams params) throws IOException {
    final File refDir = params.buildParams().directory();
    final int wordSize = params.buildParams().windowSize();
    if (params.installBlacklist() && HashBlacklist.blacklistExists(refDir, wordSize)) {
      throw new NoTalkbackSlimException("Blacklist already exists in " + refDir + " for word size " + wordSize);
    }
    final int hashBits = CreateParams.calculateHashBits(params.buildParams().sequences().mode().codeType().bits(), params.buildParams().windowSize());
    final long counterSizeBase = params.buildParams().sequences().reader().totalLength();
    final long counterSize = counterSizeBase + (long) ((params.hashMapSizeFactor() - 1.0) * counterSizeBase);
    final HashCounter sparseIndex = new HashCounter(counterSize, hashBits, params.threshold());
    final SimpleThreadPool stp = new SimpleThreadPool(params.numberThreads(), "HashToolsThread", true);
    for (int i = 0; i < params.buildParams().sequences().numberSequences(); ++i) {
      final ExactHashFunction exf = new ExactHashFunction(params.buildParams());
      final BuildParams bp = params.buildParams().subSequence(new HashingRegion(i, i + 1));
      final HashLoop subjectHashLoop = new IncrementalHashLoop(params.buildParams().stepSize(), exf, false) {
        @Override
        public void hashCall(final long hash, final int internalId, final int stepPosition) {
          //System.err.println("build hashCall hash=" + hash + " id=" + internalId);
          try {
            sparseIndex.increment(hash);
          } catch (HashCounter.TooManyCollisionsException e) {
            throw new NoTalkbackSlimException("Too many collisions in hashmap, try increasing the hashmap size factor");
          }
        }

        @Override
        public void hashCallBidirectional(long hashForward, long hashReverse, int stepPosition, int internalId) {
          throw new UnsupportedOperationException(); //"Not supported yet."
        }
      };
      executeLoop(stp, bp, subjectHashLoop);
    }
    stp.terminate();

    try (final BufferedWriter blacklistWriter = createBlacklistStream(params);
         final BufferedWriter histogramWriter = createHistogramStream(params)
    ) {
      final TreeMap<Long, Long> histMap = new TreeMap<>();
      while (sparseIndex.next()) {
        final long hash = sparseIndex.getKey();
        final long count = sparseIndex.getCount();
        if (params.makeBlacklist()) {
          if (count >= params.blacklistThreshold()) {
            blacklistWriter.write(reverseHash(hash, params.buildParams().windowSize(), params.buildParams().sequences().mode().codeType().bits()));
            blacklistWriter.write("\t");
            blacklistWriter.write(Long.toString(count));
            blacklistWriter.newLine();
          }
        }
        final long prevVal = histMap.containsKey(count) ? histMap.get(count) : 0;
        histMap.put(count, prevVal + 1);
      }
      for (Map.Entry<Long, Long> entry : histMap.entrySet()) {
        histogramWriter.append(entry.getKey().toString()).append(" ").append(entry.getValue().toString()).append(StringUtils.LS);
      }
    }
    if (params.installBlacklist()) {
      HashBlacklist.installBlacklist(params.file(BLACKLIST_FILENAME), refDir, wordSize);
    }
  }

  @SuppressWarnings("try")
  private static boolean executeLoop(SimpleThreadPool stp, BuildParams bp, HashLoop subjectHashLoop) {
    return stp.execute(() -> {
        try (BuildParams ignored = bp) {
          final OneShotTimer readTimer = new OneShotTimer("BS_read");
          subjectHashLoop.execLoop(bp.sequences(), new byte[bp.sequences().reader().length(bp.sequences().region().getStart())]);
          readTimer.stopLog();
        }
    });
  }

  private static BufferedWriter createBlacklistStream(HashDistParams params) throws IOException {
    if (params.makeBlacklist()) {
      return new BufferedWriter(new FileWriter(params.file(BLACKLIST_FILENAME)));
    }
    return null;
  }

  private static BufferedWriter createHistogramStream(HashDistParams params) throws IOException {
    return new BufferedWriter(new FileWriter(params.file("histogram.txt")));
  }

    //only nucleotides
  private static String reverseHash(long hash, int windowSize, int bitsPerThing) {
    final char[] ret = new char[windowSize];
    final int mask = (1 << bitsPerThing) - 1;
    long locHash = hash;
    for (int i = ret.length - 1; i >= 0; --i) {
      final int valuePos = 1 + ((int) locHash & mask);
      if (valuePos >= DNA.valueChars().length) {
        throw new ArrayIndexOutOfBoundsException("valuePos was: " + valuePos + " mask: " + mask + " locHash: " + locHash);
      }
      ret[i] = DNA.valueChars()[valuePos];
      locHash = locHash >>> bitsPerThing;
    }
    return new String(ret);
  }


}
