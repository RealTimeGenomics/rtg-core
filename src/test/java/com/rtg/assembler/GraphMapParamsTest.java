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

import com.rtg.assembler.graph.implementation.GraphKmerAttribute;
import com.rtg.util.IntegerOrPercentage;
import com.rtg.util.TestUtils;
import com.rtg.util.test.params.TestParams;

import junit.framework.TestCase;

/**
 */
public class GraphMapParamsTest extends TestCase {
  public void testOmnes() {
    new TestParams(GraphMapParams.class, GraphMapParams.Builder.class).check();
  }
  public void testToString() {
    final GraphMapParams params = GraphMapParams.builder().create();
    TestUtils.containsAll(params.toString()
        , "directory=" + null
        , "reads=[]"
        , "reads454=[]"
        , "readsMatePair=[]"
        , "graph=" + null
        , "wordSize=" + 0
        , "stepSize=" + 0
        , "maxMismatches=" + 0
        , "alignmentFile=" + null
        , "maxInsertSize=" + Integer.MIN_VALUE
        , "minInsertSize=" + Integer.MAX_VALUE
        , "numberThreads=" + 1
    );
  }
  public void testAssign() {
    final List<File> reads = Arrays.asList(new File("reads1"), new File("reads2"));
    final GraphMapParams params = GraphMapParams.builder()
        .directory(new File("out"))
        .graph(new GraphKmerAttribute(0))
        .reads(reads)
        .wordSize(33)
        .stepSize(11)
        .maxMismatches(new IntegerOrPercentage(4))
        .alignmentFile(new File("alignment"))
        .numberThreads(9)
        .create();
    assertEquals(new File("out"), params.directory());
    assertEquals(new File("out/foo"), params.file("foo"));
    assertNotNull(params.graph());
    assertEquals(33, params.wordSize());
    assertEquals(11, params.stepSize());
    assertEquals(9, params.numberThreads());
    assertEquals(new IntegerOrPercentage(4), params.maxMismatches());
    assertEquals(reads, params.reads());
    assertEquals(new File("alignment"), params.alignmentFile());
  }
}
