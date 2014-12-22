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

package com.rtg.assembler;

import java.io.File;
import java.util.Arrays;
import java.util.List;

import com.rtg.util.IntegerOrPercentage;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.TestUtils;
import com.rtg.util.test.params.TestParams;

import junit.framework.TestCase;

/**
 */
public class AssembleParamsTest extends TestCase {
  public void testOmnes() {
    new TestParams(AssembleParams.class, AssembleParams.Builder.class).check();
  }
  public void testToString() {
    AssembleParams params = AssembleParams.builder().create();
    TestUtils.containsAll(params.toString()
        , "directory=" + null
        , "reads=[]"
        , "reads454=[]"
        , "readsMatePair=[]"
        , "graph=" + null
        , "wordSize=" + 0
        , "stepSize=" + 0
        , "kmerSize=" + 0
        , "maxMismatches=" + 0
        , "consensusThreshold=" + 0
        , "minPathReads=" + -1
        , "minReadCount=" + -1
        , "maxMismatches=" + 0
        , "maxMismatches=" + 0
        , "alignments=" + false
        , "maxInsertSize=" + Integer.MIN_VALUE
        , "minInsertSize=" + Integer.MAX_VALUE
        , "numberThreads=" + 1
        , " region=[All inclusive]"
    );
  }
  public void testAssign() {
    List<File> reads = Arrays.asList(new File("reads1"), new File("reads2"));
    List<File> reads454 = Arrays.asList(new File("reads1a"), new File("reads2a"));
    List<File> readsMatePair = Arrays.asList(new File("reads1b"), new File("reads2b"));
    final LongRange r = new LongRange(1, 2);
    AssembleParams params = AssembleParams.builder()
        .directory(new File("out"))
        .graph(new File("in"))
        .reads(reads)
        .reads454(reads454)
        .readsMatePair(readsMatePair)
        .wordSize(33)
        .stepSize(11)
        .kmerSize(5)
        .minInsertSize(10)
        .maxInsertSize(11)
        .minPathReads(12)
        .minReadCount(13)
        .minHashFrequency(99)
        .maxMismatches(new IntegerOrPercentage(4))
        .mergeRatio(8.3)
        .alignments(true)
        .numberThreads(17)
        .region(r)
        .create();
    assertEquals(new File("out"), params.directory());
    assertEquals(new File("out/foo"), params.file("foo"));
    assertEquals(new File("in"), params.graph());
    assertEquals(33, params.wordSize());
    assertEquals(11, params.stepSize());
    assertEquals(5, params.kmerSize());
    assertEquals(10, params.minInsertSize());
    assertEquals(11, params.maxInsertSize());
    assertEquals(12, params.minPathReads());
    assertEquals(13, params.minReadCount());
    assertEquals(99, params.minHashFrequency());
    assertEquals(8.3, params.mergeRatio());
    assertEquals(17, params.numberThreads());
    assertEquals(reads454, params.reads454());
    assertEquals(readsMatePair, params.readsMatePair());
    assertEquals(new IntegerOrPercentage(4), params.maxMismatches());
    assertEquals(reads, params.reads());
    assertEquals(r, params.region());
    assertTrue(params.alignments());
  }
}
