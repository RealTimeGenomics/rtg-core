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
package com.rtg.index.hash;

import java.io.IOException;

import com.rtg.launcher.HashingRegion;
import com.rtg.launcher.ISequenceParams;
import com.rtg.mode.Frame;
import com.rtg.mode.SequenceMode;
import com.rtg.reader.SequencesReader;
import com.rtg.util.ProgramState;

/**
 * Loops over all hash codes derived from specified sequences in a reader.
 * This version increments the hash accumulator at each step.
 * It can currently only be used with <code>ExactHashFunction</code>
 */
public abstract class IncrementalHashLoop extends HashLoop {


  protected final ExactHashFunction mFunction;

  //do forward and reverse frame on same pass, assumes bidirectional frame
  protected final boolean mDualMode;

  /**
   * Constructs an instance and as a side effect makes calls to
   * <code>hashCall</code>.
   *
   * @param stepSize step size
   * @param function used to construct the hash from the codes.
   * @param dualMode whether to process forward and reverse simultaneously or not
   */
  public IncrementalHashLoop(int stepSize, final ExactHashFunction function, final boolean dualMode) {
    super(function.getWindowSize(), stepSize);
    mFunction = function;
    mDualMode = dualMode;
  }

  @Override
  public long execLoop(ISequenceParams params, final byte[] byteBuffer) throws IOException {
    return execLoop(byteBuffer, mWindowSize, mStepSize, params);
  }

  protected long execLoop(final byte[] byteBuffer, final int windowSize, final int stepSize, final ISequenceParams sequences) throws IOException {
    final HashingRegion region = sequences.region();
    final SequencesReader reader = sequences.reader();
    final long startSequence;
    final int padding = getThreadPadding();
    if (region != HashingRegion.NONE) {
      startSequence = region.getStart();
    } else {
      startSequence = 0;
    }
    final SequenceMode mode = sequences.mode();

    //number of codes to be ignored (unknown etc).
    final int firstCode = mode.codeType().firstValid();
    //the various possible frames (direction and phase for translation)
    final Frame[] frames = mode.allFrames();
    final int codeIncrement = mode.codeIncrement();
    final int jcpAdjustment = codeIncrement - windowSize * codeIncrement;
    final int step = stepSize * codeIncrement;
    //System.err.println("start=" + start + " end=" + end + " stepSize=" + stepSize + " step=" + step + " codeIncrement=" + codeIncrement);

    int internalId = (int) startSequence * frames.length;
    final long maxSequenceEver = reader.numberSequences();
    long totalLength = 0;
    for (long seq = startSequence; region.isInRange(seq) && seq < maxSequenceEver; seq++) {
      ProgramState.checkAbort();
      //System.err.println("seq=" + seq + " " +  reader.currentSequenceId() + " " + reader.getClass());
      final int currentLength = reader.length(seq);
      final int startPos = region.getReferenceStart(seq, padding);
      final int endPos = region.getReferenceEnd(seq, padding, currentLength);
      if (byteBuffer.length < endPos - startPos) { //silly programmer error
        throw new IllegalArgumentException("Allocated buffer too short. Allocated length=" + byteBuffer.length + " Required length=" + (endPos - startPos));
      }
      final int length = reader.read(seq, byteBuffer, startPos, endPos - startPos);
      totalLength += length;
      nextSeq((int) (seq - startSequence), length);
      //System.err.println("finished reading");
      //System.err.println(Arrays.toString(byteBuffer));
      for (final Frame frame : frames) {
        if (!mDualMode || frame.isForward()) {
          //System.err.println("frame=" + frame);
          next(seq, frame);
          final int phase = frame.phase();
          int limit = length;
          limit = limit - codeIncrement + 1 - phase;
          //System.err.println("limit=" + limit + " phase=" + phase + " length=" + length + " codeIncrement=" + codeIncrement);
          mFunction.reset();
          for (int j = 0; j < limit; j += codeIncrement) {
            //System.err.println("j=" + j);
            final byte b = frame.code(byteBuffer, length, j);
            final int c = b - firstCode;
            //System.err.println("b=" + b + " c=" + c);
            if (c < 0) {
              //System.err.println("reset");
              mFunction.reset();
              continue; //this code and hence window ignored
            }
            //System.err.println("j=" + j + " c=" + c);
            final long hash = mFunction.hashStep((byte) c);
            final int jcp = j + jcpAdjustment + startPos;
            //System.err.println("seq=" + seq + " frame=" + frame + " j=" + j + " jcp=" + jcp + " step=" + step);
            if (mFunction.isValid() && jcp % step == 0) {
              assert jcp >= 0;
              //System.err.println("hashCall");
              if (mDualMode) {
                hashCallBidirectional(mFunction.hash(), mFunction.hashReverse(), jcp / step, internalId);
              } else {
                hashCall(hash, internalId, jcp / step);
              }
            }
          } //window
          end();
        }
        internalId++;
      }
      //System.err.println("Finished loop");
      endSequence();

    }
    endAll();
    return totalLength;
  }

  @Override
  public void next(final long seq, final Frame frame) {
    //default do nothing
  }

  @Override
  public void end() throws IOException {
    //default do nothing
  }

  @Override
  public void endSequence() throws IOException {
    //default do nothing
  }

  @Override
  public void endAll() throws IOException {
    //default do nothing
  }
}

