/*
 * Copyright (c) 2016. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */

package com.rtg.reader;

import static org.junit.Assert.assertEquals;

import java.util.Arrays;
import java.util.List;
import java.util.function.Consumer;

import org.junit.BeforeClass;
import org.junit.Test;

import com.rtg.alignment.SingleIndelSeededEditDistance;
import com.rtg.alignment.UnidirectionalAdaptor;
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
    assertEquals(2, stats.mTotal);
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
      new UnidirectionalAdaptor(new SingleIndelSeededEditDistance(ngsParams, false, seedLength, 2, 2, maxReadLength)),
      5, 90, 4, 0, 0, false, false);
  }
}