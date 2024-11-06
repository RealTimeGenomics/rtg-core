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
package com.rtg.position.output;

import java.io.IOException;

import com.rtg.launcher.BuildParams;
import com.rtg.launcher.ISequenceParams;
import com.rtg.launcher.MockReaderParams;
import com.rtg.launcher.MockSequenceParams;
import com.rtg.mode.SequenceMode;
import com.rtg.mode.SequenceType;
import com.rtg.reader.MockSequencesReader;
import com.rtg.util.TestUtils;
import com.rtg.util.array.WrappedIntArray;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * Test GappedRegion
 */
public class GappedRegionTest extends AbstractGappedRegionTest<GappedRegion> {

  @Override
  protected GappedRegion getRegion(final int id, final int stepSize, final int wordSize, final GapScorer prob, final WrappedIntArray buildLengths) throws IOException {
    final ISequenceParams seqParams = new MockSequenceParams(new MockReaderParams(new MockSequencesReader(SequenceType.DNA)), SequenceMode.BIDIRECTIONAL);
    final BuildParams params = BuildParams.builder().windowSize(wordSize).stepSize(stepSize).sequences(seqParams).create();
    return getRegion(id, params, prob, buildLengths);
  }

  @Override
  GappedRegion getRegion(int id, BuildParams params, GapScorer prob, WrappedIntArray buildLengths) {
    return new GappedRegion(id, params, prob);
  }
  @Override
  public void setUp() {
    Diagnostic.setLogStream();
  }

  public void test() throws IOException {
    check("1[16..27]2:[3..14]^13->2", "1[0..11]0:[5..16]^5->3");
  }

  public void testNext() throws IOException {
    checkNext("2[24..35]1:[19..30]^6->1");
  }

  /** zero delta. */
  public void testMerge1() throws IOException {
    checkMerge1("2[8..35]1:[3..30]^6->3");
  }

  /** +1 delta. */
  public void testMerge2() throws IOException {
    checkMerge2a("2[8..35]1:[3..31]^6->4");
  }

  /** -1 delta. */
  public void testMerge3() throws IOException {
    final String expected = "2[8..35]1:[3..29]^6->4";
    checkMerge3(expected);
  }

  /** Overlap. */
  public void testMerge4() throws IOException {
    final String expected = "2[8..27]1:[3..22]^6->4";
    checkMerge4(expected);
  }

  /** Overlap. */
  public void testMerge6() throws IOException {
    final String expected = "2[0..23]1:[0..23]^12->4";
    checkMerge6(expected);
  }

  /** Overlap. */
  public void testMerge7() throws IOException {
    final String expected = "2[0..23]1:[0..23]^12->4";
    checkMerge7(expected);
  }

  public void testStateEnum() {
    TestUtils.testEnum(AbstractGappedRegion.State.class, "[EMPTY, INIT, WRITTEN]");
  }

}
