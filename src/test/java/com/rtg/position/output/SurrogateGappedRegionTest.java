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
import com.rtg.mode.BidirectionalFrame;
import com.rtg.mode.SequenceMode;
import com.rtg.mode.SequenceType;
import com.rtg.reader.MockSequencesReader;
import com.rtg.util.StringUtils;
import com.rtg.util.io.MemoryPrintStream;

import junit.framework.TestCase;


/**
 */
public class SurrogateGappedRegionTest extends TestCase {

  GapBuckets<GappedRegion> makeBuckets(final int stepSize) throws IOException {
    final BuildParams buildParams = GapBucketsTest.makeParams(stepSize, new int[] {30, 31, 48}, 0, SequenceMode.UNIDIRECTIONAL);
    final GapBucketsInfo bucketInfo = new GapBucketsInfo(buildParams, 1, 16, 1000);
    return new GapBuckets<>(bucketInfo);
  }

  public void testInit() throws Exception {
    final GappedDistribution distr = new GappedDistribution(2, 2, GappedDistribution.distrParams(3 * 2));
    final GapScorer prob = distr.probabilities();
    final ISequenceParams seqParams = new MockSequenceParams(new MockReaderParams(new MockSequencesReader(SequenceType.DNA)), SequenceMode.BIDIRECTIONAL);
    final BuildParams params = BuildParams.builder().windowSize(2).stepSize(2).sequences(seqParams).create();
    final GappedRegion gr = new GappedRegion(0, params, prob);
    final GapBuckets<GappedRegion> gb = makeBuckets(8);
    gr.initialize(2, 16, 3, gb, gr, false);

    final SurrogateGappedRegion sgr = new SurrogateGappedRegion(gr, 7, BidirectionalFrame.FORWARD, new BidirectionalFrame[] {BidirectionalFrame.FORWARD});

    assertTrue(sgr.scoreAllowed());
    assertEquals(2.0, sgr.score());
    assertEquals(7, sgr.queryId());
    assertEquals(sgr, sgr.initialize(gr, 7, BidirectionalFrame.REVERSE));

    final MemoryPrintStream mps = new MemoryPrintStream();
    assertTrue(sgr.write(mps.printStream(), null, null));
    assertEquals("7\tR\t4\t5\t2\tF\t17\t18" + StringUtils.LS, mps.toString());

    sgr.writeHeader(mps.printStream());
    assertTrue(mps.toString().contains("#query-id"));
  }
}
