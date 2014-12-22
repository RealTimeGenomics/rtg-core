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

import com.rtg.util.intervals.LongRange;
import com.rtg.util.TestUtils;
import com.rtg.util.test.params.TestParams;

import junit.framework.TestCase;

/**
 *         Date: 11/05/12
 *         Time: 1:30 PM
 */
public class DeBruijnParamsTest extends TestCase {
  public void testOmnes() {
    new TestParams(DeBruijnParams.class, DeBruijnParams.Builder.class).check();
  }

  public void testDefaults() {
    final DeBruijnParams def = DeBruijnParams.builder().create();
    assertEquals(0, def.kmerSize());
    assertEquals(null, def.directory());
    assertEquals(null, def.inputFiles());
    TestUtils.containsAll(def.toString()
        , "DeBruijnParams "
        , " directory=" + null
        , " inputFiles=" + null
        , " kmerSize=" + 0
        , " minHashFrequency=" + 0
        , " useStringKmers=" + false
        , " mergeRatio=0.0"
        , " region=[All inclusive]"
    );
  }
  public void testBuilder() {
    final DeBruijnParams.Builder b = DeBruijnParams.builder();
    assertEquals(b, b.kmerSize(20));
    assertEquals(b, b.inputFiles(Arrays.asList(new File("foo"), new File("bar"))));
    assertEquals(b, b.directory(new File("output")));
    assertEquals(b, b.minHashFrequency(42));
    assertEquals(b, b.useStringKmers(true));
    assertEquals(b, b.self());
    final LongRange r = new LongRange(1, 2);
    assertEquals(b, b.region(r));
    final DeBruijnParams params = b.create();
    assertEquals(20, params.kmerSize());
    assertEquals(new File("foo"), params.inputFiles().get(0));
    assertEquals(new File("output"), params.directory());
    assertEquals(42, params.minHashFrequency());
    assertEquals(r, params.region());
    assertTrue(params.useStringKmers());
  }
}
