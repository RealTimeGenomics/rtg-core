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
package com.rtg.protein;


import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.index.IndexSet;
import com.rtg.index.hash.ngs.HashFunctionFactory;
import com.rtg.index.hash.ngs.NgsHashFunction;
import com.rtg.index.hash.ngs.OutputProcessor;
import com.rtg.index.hash.ngs.ProteinNgsHashLoop;
import com.rtg.index.hash.ngs.ReadCallImplementation;
import com.rtg.index.hash.ngs.ReadEncoder;
import com.rtg.index.hash.ngs.TemplateCallImplementation;
import com.rtg.index.params.CreateParams;
import com.rtg.index.params.CreateParams.CreateParamsBuilder;
import com.rtg.launcher.HashingRegion;
import com.rtg.launcher.ISequenceParams;
import com.rtg.ngs.NgsParams;
import com.rtg.ngs.NgsTask;
import com.rtg.reader.SequencesReader;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * A couple of methods for running protein searches that I couldn't really think of anywhere better to put.
 * These functions are used instead of some corresponding NgsTask static functions.
 *         Date: 30/09/11
 *         Time: 2:44 PM
 */
@TestClass(value = {"com.rtg.protein.MapXCliTest", "com.rtg.protein.ProteinReadIndexerTest"})
public final class ProteinReadIndexer {

  /** default meta chunk size */
  public static final int DEFAULT_META_CHUNK_LENGTH = 63;
  /** default meta chunk overlap */
  public static final int DEFAULT_META_CHUNK_OVERLAP = 31;

  private static final long META_CHUNKED_KEY = -1L;

  private ProteinReadIndexer() { }
  /**
   * Construct indexes for a set of variable length protein reads.
   * @param sequence SequenceParams containing the reads to index
   * @param map of <code>bucketIds</code> to the appropriate <code>hashloop/hashfunction</code>
   * @param compressHashes should the index be compressed
   * @param numberThreads number of threads to use in indexing
   * @param buckets the read length buckets
   * @throws IOException when things have all gone horribly wrong
   * @return total number of nucleotides read.
   */
  private static long indexVariableLengths(ISequenceParams sequence, Map<Long, ReadLengthHashingState> map, boolean compressHashes, int numberThreads, SequenceLengthBuckets buckets, int metaChunkLength, int metaChunkOverlap) throws IOException {

    final long start;
    final long end;
    final HashingRegion region = sequence.region();
    final SequencesReader reader =  sequence.reader();
    if (region != HashingRegion.NONE) {
      start = region.getStart();
      end = region.getEnd();
    } else {
      start = 0;
      end = reader.numberSequences();
    }
    long totalLength = 0;
    for (int pass = 1; pass <= (compressHashes ? 2 : 1); ++pass) {
      totalLength = 0; //only use the lengths from the last pass to avoid counting input more than once.
      for (long i = start; i < end; ++i) {
        final long length = reader.length(i);
        final long bucket = buckets.bucket(length);
        if (bucket != -1) {
          final int numChunks = countMetaChunks((int) bucket, metaChunkLength, metaChunkOverlap);
          if (numChunks == 1) {
            final ReadLengthHashingState pair = map.get(bucket);
            final ProteinNgsHashLoop hashLoop = pair.getHashLoop();
            final NgsHashFunction hashFunction = pair.getHashFunction();

            final HashingRegion r = new HashingRegion(i, i + 1);
            totalLength += hashLoop.readLoop(sequence.subSequence(r), hashFunction, ReadEncoder.SINGLE_END, false);
          } else {
            final ReadLengthHashingState pair = map.get(META_CHUNKED_KEY);
            final ProteinNgsHashLoop hashLoop = pair.getHashLoop();
            final NgsHashFunction hashFunction = pair.getHashFunction();
            int readStartPosition = 0;
            for (int j = 0; j < numChunks; ++j) {
              int readEndPosition = readStartPosition + metaChunkLength * 3;
              if (readEndPosition > length) {
                readEndPosition = (int) length;
                readStartPosition = readEndPosition - metaChunkLength * 3;
              }
              final HashingRegion r = new HashingRegion(i, readStartPosition, i, readEndPosition, readStartPosition, readEndPosition);
              r.setMapxMetaChunkId(j);
              totalLength += hashLoop.readLoop(sequence.subSequence(r), hashFunction, ReadEncoder.SINGLE_END, false);
              readStartPosition = readEndPosition - metaChunkOverlap * 3;
            }
          }
        }
      }
      for (final ReadLengthHashingState indexes : map.values()) {
        final IndexSet is = indexes.getIndexes();
        is.freeze(numberThreads);
      }
    }
    return totalLength;
  }

  static ReadLengthHashingState createHashingState(int readLength, final NgsParams params, final CreateParamsBuilder indexParamsBuilder, final long count, final int frames, final OutputProcessor outProcessor, long numValues) throws IOException {
    final int windowSize = readLength / 3 - 1; //-1 to force removal of the last frame (the case where this frame is legit will just have to deal with it)
    final int maskLength = Math.min(readLength, params.mapXMetaChunkSize() * 3);
    final HashFunctionFactory hashFunctionFactory = params.maskParams().maskFactory(maskLength);
    indexParamsBuilder.size(count * frames)
    .windowBits(hashFunctionFactory.windowBits())
    .hashBits(hashFunctionFactory.hashBits());
    final CreateParams indexParams = indexParamsBuilder.create();
    final IndexSet indexes = new IndexSet(params, indexParams, hashFunctionFactory.numberWindows());
    final ProteinNgsHashLoop shl = new ProteinNgsHashLoop(windowSize, windowSize);
    Diagnostic.developerLog("index params: " + indexParams);
    final ReadCallImplementation rci = new ReadCallImplementation(indexes);
    final TemplateCallImplementation tci = new TemplateCallImplementation(params, numValues, indexes, outProcessor);
    final NgsHashFunction hf = hashFunctionFactory.create(rci, tci);
    final long numberReads = params.buildFirstParams().numberSequences();
    hf.setReadSequences(numberReads);
    return new ReadLengthHashingState(indexes, shl, hf, tci);
  }

  /**
   * builds an index then runs a search
   * @param params search parameters
   * @param outProcessor output destination
   * @param indexParamsBuilder index creation parameter builder.
   * @param buckets the sequence length divisions
   * @param numValues range of different values that can be stored against a hash
   * @throws IOException when output doesn't want to play
   * @return the count of the number of residues read.
   */
  public static long indexThenSearchProteinReads(final NgsParams params, final OutputProcessor outProcessor, final CreateParams.CreateParamsBuilder indexParamsBuilder, SequenceLengthBuckets buckets, long numValues) throws IOException {
    assert !params.paired();
    final Map<Long, ReadLengthHashingState> lengthFunctions = new HashMap<>();
    final int frames = params.buildFirstParams().mode().numberFrames();

    long longReadCount = 0;
    for (final Long l : buckets.getBuckets()) {
      final int numChunks = countMetaChunks(l.intValue(), params.mapXMetaChunkSize(), params.mapXMetaChunkOverlap());
      if (numChunks > 1) {
        longReadCount += numChunks * buckets.getCount(l);
      } else {
        final long count = buckets.getCount(l);
        final ReadLengthHashingState rlhs = createHashingState(l.intValue(), params, indexParamsBuilder, count, frames, outProcessor, numValues);
        lengthFunctions.put(l, rlhs);
      }
    }
    if (longReadCount > 0) {
      final int virtualReadLength = params.mapXMetaChunkSize() * 3;
      final ReadLengthHashingState rlhs = createHashingState(virtualReadLength, params, indexParamsBuilder, longReadCount, frames, outProcessor, numValues);
      lengthFunctions.put(META_CHUNKED_KEY, rlhs);
    }
    final long totalLength = ProteinReadIndexer.indexVariableLengths(params.buildFirstParams(), lengthFunctions, params.compressHashes(), params.numberThreads(), buckets, params.mapXMetaChunkSize(), params.mapXMetaChunkOverlap());
    int num = 0;
    final int total = lengthFunctions.size();
    for (final Map.Entry<Long, ReadLengthHashingState> entry : lengthFunctions.entrySet()) {
      //Diagnostic.developerLog("Searching for read length: " + l + " total length groups: " + buckets.getBuckets().size());
      Diagnostic.progress("ProteinSearchBucket: " + ++num + "/" + total);
      final ReadLengthHashingState hashFunctions = entry.getValue();
      NgsTask.search(params, hashFunctions.getHashLoop(), hashFunctions.getIndexes(), hashFunctions.getTemplateCallImplementation(), hashFunctions.getHashFunction());
    }
    return totalLength;
  }

  /**
   * Returns start positions (protein space) for each chunk given DNA read length and overlap size
   * @param length the length of the read in DNA space
   * @param metaChunkLength size of meta chunks
   * @param metaChunkOverlap  overlap size in protein space
   * @return start position for each chunk in protein space
   */
  static int countMetaChunks(int length, int metaChunkLength, int metaChunkOverlap) {
    final int protlength = (length / 3) - 1;
    //63 = max protein length
    return (protlength - metaChunkOverlap - 1) / (metaChunkLength - metaChunkOverlap) + 1;
  }

  static int chunkToPosition(int chunkId, int plen, int metaChunkLength, int metaChunkOverlap) {
    if (chunkId == 0) {
      return 0;
    }
    final int ret = chunkId * (metaChunkLength - metaChunkOverlap);
    if (ret > plen - metaChunkLength) {
      return plen - metaChunkLength;
    }
    return ret;
  }
}
