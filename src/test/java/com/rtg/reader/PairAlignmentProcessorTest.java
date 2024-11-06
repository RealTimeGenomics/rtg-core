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

package com.rtg.reader;

import static org.junit.Assert.assertEquals;

import java.util.Arrays;
import java.util.List;
import java.util.function.Consumer;

import org.junit.BeforeClass;
import org.junit.Test;

import com.rtg.alignment.SingleIndelSeededEditDistance;
import com.rtg.ngs.NgsParams;
import com.rtg.ngs.NgsParamsBuilder;
import com.rtg.util.diagnostic.Diagnostic;

public class PairAlignmentProcessorTest {
  @BeforeClass
  public static void setUp() {
    Diagnostic.setLogStream();
  }


  @Test
  public void testAlignmentProcessor() {
    final StringBuilder r1Sb = new StringBuilder();
    final StringBuilder r2Sb = new StringBuilder();

    final Consumer<List<FastqPair>> consumer = list -> list.forEach(item -> {
      r1Sb.append(item.r1().toFasta());
      r2Sb.append(item.r2().toFasta());
    });
    final Batch<FastqPair> batch = new Batch<>(0, Arrays.asList(
      getFastqPair("r1", "GACGACGTTTGTT", "AACAAACGTCGTC"),
      getFastqPair("r2", "GACGACGAAAGTTGG", "AACTTTCGTCGTCCC")

    ));
    final PairAlignmentStats stats = new PairAlignmentStats();
    final PairAlignmentProcessor pairAlignmentProcessor = new PairAlignmentProcessor(stats, new BatchReorderingWriter<>(consumer), batch, getPairAligner());
    pairAlignmentProcessor.run();
    final String expectedR1 = fasta("r1", "GACGACGTTTGTT")
      + fasta("r2", "GACGACGAAAGTT");
    final String expectedR2 = fasta("r1", "AACAAACGT")
      + fasta("r2", "AACTTTCGT");
    assertEquals(expectedR1, r1Sb.toString());
    assertEquals(expectedR2, r2Sb.toString());
    assertEquals(1, stats.mR2ReadIntoR1Probe);
    assertEquals(1, stats.mR1ReadThrough);
    assertEquals(1, stats.mR2ReadThrough);
    assertEquals(2, stats.mTotalInput);
  }
  private String fasta(String name, String sequence) {
    return ">" + name + "\n" + sequence + "\n";
  }

  private FastqPair getFastqPair(String name, String  r1, String r2) {
    return new FastqPair(FastqSequenceTest.getFastq(name, r1), FastqSequenceTest.getFastq(name, r2));
  }

  private PairAligner getPairAligner() {
    final int maxReadLength = 300;
    final int seedLength = 5;

    final NgsParams ngsParams = new NgsParamsBuilder().create();
    return new PairAligner(
      new SingleIndelSeededEditDistance(ngsParams, false, seedLength, 2, 2, maxReadLength),
      5, 90, 4, 0, 0, false, false, PairAligner.MismatchType.NONE, false);
  }
}
