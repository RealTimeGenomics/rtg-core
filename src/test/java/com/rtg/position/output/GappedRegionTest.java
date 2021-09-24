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
