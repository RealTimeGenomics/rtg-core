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
public class HashDist extends ParamsTask<HashDistParams, NullStatistics> {

  private static final String BLACKLIST_FILENAME = "blacklist";

  /**
   * @param params parameters
   * @param reportStream unused
   * @param stats unused
   * @param usageMetric unused
   */
  public HashDist(HashDistParams params, OutputStream reportStream, NullStatistics stats, UsageMetric usageMetric) {
    super(params, reportStream, stats, usageMetric);
  }

  @Override
  protected void exec() throws IOException {
    hashCount(mParams);
  }

  static void hashCount(final HashDistParams params) throws IOException {
    final File refDir = params.build().directory();
    final int wordSize = params.build().windowSize();
    if (params.installBlacklist() && HashBlacklist.blacklistExists(refDir, wordSize)) {
      throw new NoTalkbackSlimException("Blacklist already exists in " + refDir + " for word size " + wordSize);
    }
    final int hashBits = CreateParams.calculateHashBits(params.build().sequences().mode().codeType().bits(), params.build().windowSize());
    final long counterSizeBase = params.build().sequences().reader().totalLength();
    final long counterSize = counterSizeBase + (long) ((params.hashMapSizeFactor() - 1.0) * counterSizeBase);
    final HashCounter sparseIndex = new HashCounter(counterSize, hashBits, params.threshold());
    final SimpleThreadPool stp = new SimpleThreadPool(params.numberThreads(), "HashToolsThread", true);
    for (int i = 0; i < params.build().sequences().numberSequences(); i++) {
      final ExactHashFunction exf = new ExactHashFunction(params.build());
      final BuildParams bp = params.build().subSequence(new HashingRegion(i, i + 1));
      final HashLoop subjectHashLoop = new IncrementalHashLoop(params.build().stepSize(), exf, false) {
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
            blacklistWriter.write(reverseHash(hash, params.build().windowSize(), params.build().sequences().mode().codeType().bits()));
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
    for (int i = ret.length - 1; i >= 0; i--) {
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
