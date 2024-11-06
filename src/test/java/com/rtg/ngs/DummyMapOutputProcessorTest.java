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
package com.rtg.ngs;

import java.util.Arrays;

import com.rtg.launcher.HashingRegion;

import junit.framework.TestCase;

/**
 */
public class DummyMapOutputProcessorTest extends TestCase {

  public void testSingleResult() {
    //region.NONE (all encompassing region)
    final MatchResult results = new MatchResult(0);
    results.addMatchResult(1, 50, 0, false);
    HashingRegion[] regions = {HashingRegion.NONE};
    AbstractMapOutputProcessor.ChunkPair[] chunks;
    AbstractMapOutputProcessor.ChunkPair[] expected;
    chunks = AbstractMapOutputProcessor.findChunkBoundaries(regions, results);
    expected = new AbstractMapOutputProcessor.ChunkPair[] {new AbstractMapOutputProcessor.ChunkPair(0, 1)};
    assertTrue(Arrays.equals(expected, chunks));

    //Single region including result
    regions = new HashingRegion[] {new HashingRegion(1, 40, 1, 100, 20, 150)};
    chunks = AbstractMapOutputProcessor.findChunkBoundaries(regions, results);
    assertTrue(Arrays.equals(expected, chunks));

    //two regions, only second including result
    regions = new HashingRegion[] {new HashingRegion(0, 0, 0, 500, 0, 700), new HashingRegion(1, 40, 1, 100, 20, 150)};
    chunks = AbstractMapOutputProcessor.findChunkBoundaries(regions, results);
    expected = new AbstractMapOutputProcessor.ChunkPair[] {new AbstractMapOutputProcessor.ChunkPair(0, 0), new AbstractMapOutputProcessor.ChunkPair(0, 1)};
    assertTrue("Expected:\n" + Arrays.toString(expected) + "\nActual:\n" + Arrays.toString(chunks), Arrays.equals(expected, chunks));

    //two regions, both including result
    regions = new HashingRegion[] {new HashingRegion(1, 20, 1, 80, 10, 200), new HashingRegion(1, 40, 1, 100, 20, 150)};
    chunks = AbstractMapOutputProcessor.findChunkBoundaries(regions, results);
    expected = new AbstractMapOutputProcessor.ChunkPair[] {new AbstractMapOutputProcessor.ChunkPair(0, 1), new AbstractMapOutputProcessor.ChunkPair(0, 1)};
    assertTrue("Expected:\n" + Arrays.toString(expected) + "\nActual:\n" + Arrays.toString(chunks), Arrays.equals(expected, chunks));

    //two regions, only first including result
    regions = new HashingRegion[] {new HashingRegion(1, 20, 1, 80, 10, 200), new HashingRegion(2, 40, 2, 100, 20, 150)};
    chunks = AbstractMapOutputProcessor.findChunkBoundaries(regions, results);
    expected = new AbstractMapOutputProcessor.ChunkPair[] {new AbstractMapOutputProcessor.ChunkPair(0, 1), new AbstractMapOutputProcessor.ChunkPair(1, 1)};
    assertTrue("Expected:\n" + Arrays.toString(expected) + "\nActual:\n" + Arrays.toString(chunks), Arrays.equals(expected, chunks));

    //three regions, only third including result
    regions = new HashingRegion[] {new HashingRegion(0, 0, 0, 500, 0, 700), new HashingRegion(0, 500, 0, 600, 30, 1000), new HashingRegion(1, 40, 1, 100, 20, 150)};
    chunks = AbstractMapOutputProcessor.findChunkBoundaries(regions, results);
    expected = new AbstractMapOutputProcessor.ChunkPair[] {new AbstractMapOutputProcessor.ChunkPair(0, 0), new AbstractMapOutputProcessor.ChunkPair(0, 0), new AbstractMapOutputProcessor.ChunkPair(0, 1)};
    assertTrue("Expected:\n" + Arrays.toString(expected) + "\nActual:\n" + Arrays.toString(chunks), Arrays.equals(expected, chunks));
  }

  public void testMultiResult() {
    //region.NONE (all encompassing region)
    final MatchResult results = new MatchResult(0);
    results.addMatchResult(1, 50, 0, false);
    results.addMatchResult(2, 100, 0, false);
    results.addMatchResult(3, 500, 0, false);
    HashingRegion[] regions = {HashingRegion.NONE};
    AbstractMapOutputProcessor.ChunkPair[] chunks;
    AbstractMapOutputProcessor.ChunkPair[] expected;
    chunks = AbstractMapOutputProcessor.findChunkBoundaries(regions, results);
    expected = new AbstractMapOutputProcessor.ChunkPair[] {new AbstractMapOutputProcessor.ChunkPair(0, 3)};
    assertTrue(Arrays.equals(expected, chunks));

    //Single region including second result
    regions = new HashingRegion[] {new HashingRegion(2, 40, 2, 100, 20, 150)};
    chunks = AbstractMapOutputProcessor.findChunkBoundaries(regions, results);
    expected = new AbstractMapOutputProcessor.ChunkPair[] {new AbstractMapOutputProcessor.ChunkPair(1, 2)};
    assertTrue("Expected:\n" + Arrays.toString(expected) + "\nActual:\n" + Arrays.toString(chunks), Arrays.equals(expected, chunks));

    //three regions, second and third including overlapping results
    regions = new HashingRegion[] {new HashingRegion(0, 0, 0, 500, 0, 700), new HashingRegion(1, 40, 1, 100, 20, 150), new HashingRegion(1, 60, 2, 150, 40, 200)};
    chunks = AbstractMapOutputProcessor.findChunkBoundaries(regions, results);
    expected = new AbstractMapOutputProcessor.ChunkPair[] {new AbstractMapOutputProcessor.ChunkPair(0, 0), new AbstractMapOutputProcessor.ChunkPair(0, 1), new AbstractMapOutputProcessor.ChunkPair(0, 2)};
    assertTrue("Expected:\n" + Arrays.toString(expected) + "\nActual:\n" + Arrays.toString(chunks), Arrays.equals(expected, chunks));

    //three regions, first and second including overlapping results
    regions = new HashingRegion[] {new HashingRegion(1, 40, 1, 100, 20, 150), new HashingRegion(1, 60, 2, 150, 40, 200), new HashingRegion(4, 0, 4, 100, 0, 1000)};
    chunks = AbstractMapOutputProcessor.findChunkBoundaries(regions, results);
    expected = new AbstractMapOutputProcessor.ChunkPair[] {new AbstractMapOutputProcessor.ChunkPair(0, 1), new AbstractMapOutputProcessor.ChunkPair(0, 2), new AbstractMapOutputProcessor.ChunkPair(3, 3)};
    assertTrue("Expected:\n" + Arrays.toString(expected) + "\nActual:\n" + Arrays.toString(chunks), Arrays.equals(expected, chunks));

    //three regions, first has result, second has none, third has another result
    regions = new HashingRegion[] {new HashingRegion(2, 90, 2, 500, 80, 600), new HashingRegion(2, 1000, 2, 2000, 900, 2500), new HashingRegion(3, 490, 3, 1000, 450, 1200)};
    chunks = AbstractMapOutputProcessor.findChunkBoundaries(regions, results);
    expected = new AbstractMapOutputProcessor.ChunkPair[] {new AbstractMapOutputProcessor.ChunkPair(1, 2), new AbstractMapOutputProcessor.ChunkPair(2, 2), new AbstractMapOutputProcessor.ChunkPair(2, 3)};
    assertTrue("Expected:\n" + Arrays.toString(expected) + "\nActual:\n" + Arrays.toString(chunks), Arrays.equals(expected, chunks));

    //four regions, first has result, second and third has none, fourth has another result
    regions = new HashingRegion[] {new HashingRegion(2, 90, 2, 500, 80, 600), new HashingRegion(2, 1000, 2, 2000, 900, 2500), new HashingRegion(2, 20000, 2, 20000, 19000, 25000), new HashingRegion(3, 490, 3, 1000, 450, 1200)};
    chunks = AbstractMapOutputProcessor.findChunkBoundaries(regions, results);
    expected = new AbstractMapOutputProcessor.ChunkPair[] {new AbstractMapOutputProcessor.ChunkPair(1, 2), new AbstractMapOutputProcessor.ChunkPair(2, 2), new AbstractMapOutputProcessor.ChunkPair(2, 2), new AbstractMapOutputProcessor.ChunkPair(2, 3)};
    assertTrue("Expected:\n" + Arrays.toString(expected) + "\nActual:\n" + Arrays.toString(chunks), Arrays.equals(expected, chunks));
  }

//  public void test() throws IOException {
//    try (TestDirectory tdir = new TestDirectory()) {
//      final File dir = new File("/rtgshare/data/human/ref/1000g_v37_phase2/sdf/");
//      SequencesReader dsr = SequencesReaderFactory.createDefaultSequencesReader(dir);
//      List<HashingRegion> regions = Arrays.asList(HashingRegion.splitWorkload(dsr, Sex.MALE, 0, dsr.numberSequences(), 16 * HashingRegion.DEFAULT_THREAD_MULTIPLIER, HashingRegion.DEFAULT_MIN_CHUNK_SIZE, 1202));
//      ISequenceParams searchParam = SequenceParams.builder().directory(dir).useMemReader(true).create();
//
//      Interval[] loci = AbstractMulticoreFilterConcat.groupRegions(regions.size(), AbstractMulticoreFilterConcat.numberIntermediateFiles(regions.size(), 16));
//      UnmappedSamAlignmentWriter[] writers = new UnmappedSamAlignmentWriter[loci.length];
//      for (int i = 0; i < loci.length; ++i) {
//        writers[i] = new UnmappedSamAlignmentWriter(tdir, new SAMFileHeader());
//      }
//      final Map<Long,RangeList<UnmappedSamAlignmentWriter>> longRangeSearchMap = AbstractMapOutputProcessor.getUnmappedWriterReferenceLookup(searchParam, regions, loci, writers, new UnmappedSamAlignmentWriter(tdir, new SAMFileHeader()));
//      RangeList<UnmappedSamAlignmentWriter> moo = longRangeSearchMap.get(23L);
//      System.err.println("foo");
//    }
//  }
}
