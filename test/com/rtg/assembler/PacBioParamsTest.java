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

import com.rtg.util.TestUtils;
import com.rtg.util.test.params.TestParams;

import junit.framework.TestCase;

/**
 */
public class PacBioParamsTest extends TestCase {
  public void testOmnes() {
    new TestParams(PacBioParams.class, PacBioParams.Builder.class).check();
  }

  public void testToString() {
    final PacBioParams params = PacBioParams.builder().create();
    TestUtils.containsAll(params.toString()
        , " directory=" + null
        , " reads=" + null
        , " graph=" + null
        , " trimGraph=" + false
    );
  }

  public void testAssign() {
    final List<File> reads = Arrays.asList(new File("foo"));
    final File graph = new File("graph");
    final File out = new File("out");
    final PacBioParams pacBioParams = PacBioParams.builder().directory(out).reads(reads).graph(graph).trimGraph(true).create();
    assertEquals(graph, pacBioParams.graph());
    assertEquals(reads, pacBioParams.reads());
    assertEquals(out, pacBioParams.directory());
    assertTrue(pacBioParams.trimGraph());
    assertEquals(new File(out, "bar"), pacBioParams.file("bar"));
  }
}
