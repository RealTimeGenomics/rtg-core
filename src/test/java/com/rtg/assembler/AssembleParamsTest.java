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
    final AssembleParams params = AssembleParams.builder().create();
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
    final List<File> reads = Arrays.asList(new File("reads1"), new File("reads2"));
    final List<File> reads454 = Arrays.asList(new File("reads1a"), new File("reads2a"));
    final List<File> readsMatePair = Arrays.asList(new File("reads1b"), new File("reads2b"));
    final LongRange r = new LongRange(1, 2);
    final AssembleParams params = AssembleParams.builder()
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
