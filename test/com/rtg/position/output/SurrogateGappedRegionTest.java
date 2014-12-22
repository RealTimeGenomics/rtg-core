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
    final ISequenceParams seqParams = new MockSequenceParams(new MockReaderParams(new MockSequencesReader(SequenceType.DNA), SequenceMode.BIDIRECTIONAL));
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
