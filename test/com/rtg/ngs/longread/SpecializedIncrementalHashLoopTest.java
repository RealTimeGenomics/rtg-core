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
package com.rtg.ngs.longread;

import static com.rtg.util.StringUtils.LS;

import java.io.File;

import com.rtg.index.Index;
import com.rtg.index.hash.ExactHashFunction;
import com.rtg.launcher.HashingRegion;
import com.rtg.launcher.ISequenceParams;
import com.rtg.launcher.MockReaderParams;
import com.rtg.launcher.MockSequenceParams;
import com.rtg.launcher.ReaderParams;
import com.rtg.mode.BidirectionalFrame;
import com.rtg.mode.DnaUtils;
import com.rtg.mode.Frame;
import com.rtg.mode.SequenceMode;
import com.rtg.ngs.NgsLongTest;
import com.rtg.position.FinderPositionOutput;
import com.rtg.position.MockPositionOutput;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 */
public class SpecializedIncrementalHashLoopTest extends NgsLongTest {

  public static Test suite() {
    final TestSuite suite = new TestSuite();

    suite.addTestSuite(SpecializedIncrementalHashLoopTest.class);
    return suite;
  }

  //private static final String READ_1 = ">read_l/1" + LS
  //+ "TGCCAC" + LS
  //+ ">read_2/2" + LS + "TGCTGT" + LS;
  //12345678901234567890123456789012345678901234567890123456789012

  private static final String TEMPSTR = "CAGGCAACTGCCACCTTGGTTTTTTGCCCCCCTG";
  private static final byte[] TEMPSTRBYTES = DnaUtils.encodeString(TEMPSTR);
  private static final String TEMP_1 = ">template" + LS
  + "CAGGCAACTGCCACCTT" + LS
  + "GGTTTTTTGCCCCCCTG" + LS
  + ">template2" + LS
  + TEMPSTR + LS;

  private static class MockHashFunction1 extends ExactHashFunction {
    int mTmppos = 1;

    MockHashFunction1(final int windowSize, final int bits) {
      super(windowSize, bits);
    }

    @Override
    public long hashStep(final byte code) {
      assertEquals(TEMPSTRBYTES[mTmppos++] - 1, code);
      return super.hashStep(code);
    }
  }

  public void testRegion() throws Exception {
    final File tmpDir = FileUtils.createTempDir("longreadtask", "tmpdir");
    try {

      final File template = new File(tmpDir, "template");
      ReaderTestUtils.getReaderDNA(TEMP_1, template, null);
      final SequencesReader tmplreader = SequencesReaderFactory.createMemorySequencesReader(template, true, LongRange.NONE);

      final FinderPositionOutput fpo = new FinderPositionOutput(null, new MockPositionOutput());

      final int[] pos = {1};
      final SpecializedIncrementalHashLoop sihl = new SpecializedIncrementalHashLoop(2, new MockHashFunction1(2, 2), fpo, (FinderPositionOutput) null, (Index) null, true) {

        @Override
        public void hashCallBidirectional(final long hashFoward, final long hashReverse, final int stepPosition, final int internalId) {
          assertEquals(pos[0]++, stepPosition);
          assertEquals(0, internalId);
        }

        @Override
        public void next(final long seq, final Frame frame) {
          assertEquals(1, seq);
          mOutput.nextQuery(BidirectionalFrame.FORWARD, (int) seq);
        }
      };

      final ReaderParams rp = new MockReaderParams(tmplreader, SequenceMode.BIDIRECTIONAL);

      final ISequenceParams seqparam = new MockSequenceParams(rp) {
        @Override
        public HashingRegion region() {
          return new HashingRegion(1, 1, 2, -1, 1, -1);
        }
      };
      try {
        sihl.execLoop(new byte[1], 2, 1, seqparam);
        fail();
      } catch (final IllegalArgumentException iae) {
        assertEquals("Allocated buffer too short. Allocated length=1 Required length=33", iae.getMessage());
      }

      sihl.execLoop(new byte[33], 2, 1, seqparam);
      assertEquals(33, pos[0]);

    } finally {
      FileHelper.deleteAll(tmpDir);
    }
  }
  private static class MockHashFunction2 extends ExactHashFunction {
    MockHashFunction2(final int windowSize, final int bits) {
      super(windowSize, bits);
    }
  }

  public void testNoRegion() throws Exception {
    final File tmpDir = FileUtils.createTempDir("longreadtask", "tmpdir");
    try {

      final File template = new File(tmpDir, "template");
      ReaderTestUtils.getReaderDNA(TEMP_1, template, null);
      final SequencesReader tmplreader = SequencesReaderFactory.createMemorySequencesReader(template, true, LongRange.NONE);

      final FinderPositionOutput fpo = new FinderPositionOutput(null, new MockPositionOutput());



      final SpecializedIncrementalHashLoop sihl = new SpecializedIncrementalHashLoop(2, new MockHashFunction2(2, 2), fpo, (FinderPositionOutput) null, (Index) null, true) {
        boolean mSeenFirst = false;

        @Override
        public void hashCallBidirectional(final long hashFoward, final long hashReverse, final int stepPosition, final int internalId) {
          assertEquals(0, internalId);
        }

        @Override
        public void next(final long seq, final Frame frame) {
          if (seq == 0) {
            mSeenFirst = true;
          } else {
            assertEquals(1, seq);
            assertTrue(mSeenFirst);
          }
          mOutput.nextQuery(BidirectionalFrame.FORWARD, (int) seq);
        }
      };

      final ReaderParams rp = new MockReaderParams(tmplreader, SequenceMode.BIDIRECTIONAL);

      final ISequenceParams seqparam = new MockSequenceParamsImpl(rp);
      sihl.execLoop(new byte[1024], 2, 1, seqparam);
    } finally {
      FileHelper.deleteAll(tmpDir);
    }
  }

  private static class MockSequenceParamsImpl extends MockSequenceParams {

    public MockSequenceParamsImpl(ReaderParams readerParams) {
      super(readerParams);
    }

    @Override
    public HashingRegion region() {
      return HashingRegion.NONE;
    }
  }
}
