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

import java.io.IOException;

import com.rtg.index.Index;
import com.rtg.index.hash.ExactHashFunction;
import com.rtg.launcher.ISequenceParams;
import com.rtg.launcher.HashingRegion;
import com.rtg.mode.BidirectionalFrame;
import com.rtg.mode.Frame;
import com.rtg.mode.SequenceMode;
import com.rtg.position.FinderPositionOutput;
import com.rtg.position.SearchIncrementalHashLoop;
import com.rtg.reader.SequencesReader;
import com.rtg.util.ProgramState;

/**
 * Overrides <code>execLoop</code> for the special case of long reads.
 * This is intended to optimize the case where there are few hits
 * which is important for some industry standard benchmarks.
 */
public class SpecializedIncrementalHashLoop extends SearchIncrementalHashLoop {

  /**
   * @param stepSize step size
   * @param function hash function.
   * @param outputVars variables for output parameters
   * @param outputVarsReverse variables for output parameters (for reverse case)
   * @param index to be updated.
   * @param dualMode whether to process forward and reverse simultaneously or not
   */
  public SpecializedIncrementalHashLoop(int stepSize, final ExactHashFunction function,
                                        final FinderPositionOutput outputVars, final FinderPositionOutput outputVarsReverse, final Index index, final boolean dualMode) {
    super(stepSize, function, outputVars, outputVarsReverse, index, dualMode);
  }

  @Override
  protected long execLoop(final byte[] byteBuffer, final int windowSize, final int stepSize, final ISequenceParams sequences) throws IOException {
    assert stepSize == 1;
    assert mDualMode;
    final HashingRegion region = sequences.region();
    final SequencesReader reader = sequences.reader();
    final long startSequence;
    final long endSequence;
    final int padding = getThreadPadding();
    if (region != HashingRegion.NONE) {
      startSequence = region.getStart();
      endSequence = region.getEnd();
    } else {
      startSequence = 0;
      endSequence = reader.numberSequences();
    }
    final SequenceMode mode = sequences.mode();
    assert mode == SequenceMode.BIDIRECTIONAL;
    assert mode.allFrames().length == 2;
    assert mode.allFrames()[0] == BidirectionalFrame.FORWARD;
    assert mode.codeIncrement() == 1;
    assert mode.codeType().firstValid() == 1;

    final int jcpAdjustment = 1 - windowSize;
    //System.err.println("start=" + start + " end=" + end + " stepSize=" + stepSize + " step=" + step + " codeIncrement=" + codeIncrement);

    final long maxSequenceEver = reader.numberSequences();
    if (startSequence <= endSequence && startSequence < maxSequenceEver) {
      reader.seek(startSequence);
    }
    long totalLength = 0;
    for (long seq = startSequence; region.isInRange(seq) && seq < maxSequenceEver; seq++) {
      ProgramState.checkAbort();
      //System.err.println("seq=" + seq + " " +  reader.currentSequenceId() + " " + reader.getClass());
      final int startPos = region.getReferenceStart(seq, padding);
      final int jstart = jcpAdjustment + startPos;
      final int endPos = region.getReferenceEnd(seq, padding, reader.currentLength());
      if (byteBuffer.length < endPos - startPos) { //silly programmer error
        throw new IllegalArgumentException("Allocated buffer too short. Allocated length=" + byteBuffer.length + " Required length=" + (endPos - startPos));
      }
      final int length = reader.readCurrent(byteBuffer, startPos, endPos - startPos);
      totalLength += length;
      nextSeq((int) seq, reader.currentLength(), length, byteBuffer);
      //System.err.println("finished reading");
      //System.err.println(Arrays.toString(byteBuffer));
      //System.err.println("frame=" + frame);
      next(seq, BidirectionalFrame.FORWARD);
      //System.err.println("limit=" + limit + " phase=" + phase + " length=" + length + " codeIncrement=" + codeIncrement);
      mFunction.reset();
      for (int j = 0; j < length; j++) {
        //System.err.println("j=" + j);
        final byte b = BidirectionalFrame.FORWARD.code(byteBuffer, length, j);
        final int c = b - 1;
        //System.err.println("b=" + b + " c=" + c);
        if (c < 0) {
          //System.err.println("reset");
          mFunction.reset();
          continue; //this code and hence window ignored
        }
        //System.err.println("c=" + c);
        mFunction.hashStep((byte) c);
        //System.err.println("seq=" + seq + " frame=" + frame + " j=" + j + " jcp=" + jcp + " step=" + step);
        if (mFunction.isValid()) {
          //System.err.println("hashCall");
          hashCallBidirectional(mFunction.hash(), mFunction.hashReverse(), j + jstart, 0/*Not needed*/);
        }
      } //window
      end();
      //System.err.println("Finished loop");
      reader.nextSequence();
      endSequence();

    }
    endAll();
    return totalLength;
  }

  private void nextSeq(int seqId, final int length, final int usedLength, final byte[] read) {
    mOutput.nextSequence(seqId, length, usedLength, read);
    if (mOutputReverse != null) {
      mOutputReverse.nextSequence(seqId, length, usedLength, read);
    }
  }


  @Override
  public void next(final long seq, final Frame frame) {
    assert frame == BidirectionalFrame.FORWARD;
    mOutput.nextQuery(BidirectionalFrame.FORWARD, (int) seq);
    if (mOutputReverse != null) {
      mOutputReverse.nextQuery(BidirectionalFrame.REVERSE, (int) seq);
    }
  }
}
