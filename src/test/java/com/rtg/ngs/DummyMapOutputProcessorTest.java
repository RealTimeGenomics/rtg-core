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
