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

import com.rtg.util.intervals.LongRange;
import com.rtg.util.TestUtils;
import com.rtg.util.test.params.TestParams;

import junit.framework.TestCase;

/**
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
