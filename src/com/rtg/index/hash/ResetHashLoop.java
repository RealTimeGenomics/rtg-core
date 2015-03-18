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
 * This version resets the hash accumulator at the start of each window.
 */
public abstract class ResetHashLoop extends HashLoop {


  protected final HashFunction mFunction;

  private final boolean mDualMode;

  /**
   * Constructs an instance and as a side effect makes calls to
   * <code>hashCall</code>.
   * @param windowSize number of codes to be included in a window used to make a hash.
   * @param stepSize step size
   * @param function used to construct the hash from the codes.
   * @param dualMode whether to process forward and reverse simultaneously or not
   */
  public ResetHashLoop(final int windowSize, int stepSize, final HashFunction function, final boolean dualMode) {
    super(windowSize, stepSize);
    mFunction = function;
    mDualMode = dualMode;
  }


  @Override
  public long execLoop(ISequenceParams params, byte[] byteBuffer) throws IOException {
    return execLoop(byteBuffer, mStepSize, params);
  }

  private long execLoop(final byte[] byteBuffer, final int stepSize,
      final ISequenceParams sequences) throws IOException {
    final HashingRegion region = sequences.region();
    final SequencesReader reader = sequences.reader();
    final long startSequence;
    final int padding = getThreadPadding();
    if (region == HashingRegion.NONE) {
      startSequence = 0;
    } else {
      startSequence = region.getStart();
    }
    final SequenceMode mode = sequences.mode();
    //assert reader.integrity();
    //number of codes to be ignored (unknown etc).
    final int firstCode = mode.codeType().firstValid();
    //the various possible frames (direction and phase for translation)
    final Frame[] frames = mode.allFrames();
    final int codeIncrement = mode.codeIncrement();
    final int step = stepSize * codeIncrement;

    int internalId = (int) startSequence * frames.length;
    final long maxSequenceEver = reader.numberSequences();
    long totalLength = 0;
    //System.err.println("start=" + start + " end=" + end + " stepSize=" + stepSize + " step=" + step + " codeIncrement=" + codeIncrement);
    for (long seq = startSequence; seq < maxSequenceEver && region.isInRange(seq); seq++) {
      //System.err.println("seq=" + seq);
      final int currentLength = reader.length(seq);
      final int startPos = region.getReferenceStart(seq, padding);
      final int endPos = region.getReferenceEnd(seq, padding, currentLength);
      if (byteBuffer.length < endPos - startPos) {
        throw new IllegalArgumentException("Allocated buffer too short. Allocated length=" + byteBuffer.length + " Required length=" + (endPos - startPos));
      }
      final int length = reader.read(seq, byteBuffer, startPos, endPos - startPos);
      totalLength += length;
      nextSeq((int) (seq - startSequence), currentLength);
      //System.err.println(Arrays.toString(byteBuffer));
      final int limitOffset = mWindowSize * codeIncrement;
      final int limit0 = length - limitOffset;
      for (final Frame frame : frames) {
        if (!(mDualMode && !frame.isForward())) {
          //System.err.println("frame=" + frame);
          next(seq, frame);
          final int phase = frame.phase();
          final int limit = limit0 - phase;
          boolean cont = false;
          for (int j = 0, s = 0; j <= limit; j += step, s++) {
            ProgramState.checkAbort();
            //System.err.println("j=" + j);
            mFunction.reset();
            long hash = 0L;
            for (int kk = 0, k = j; kk < mWindowSize; kk++, k += codeIncrement) {
              //System.err.println("kk=" + kk);
              final byte b = frame.code(byteBuffer, length, k);
              final int c = b - firstCode;
              //System.err.println("b=" + b + " c=" + c);
              if (c < 0) {
                cont = true;
                break;
                //continue window; //this code and hence window ignored
              }
              //System.err.println("c=" + c);
              hash = mFunction.hashStep((byte) c);
            }
            if (cont) {
              cont = false;
              continue;
            }
            if (mDualMode) {
              final long reverseHash = reverseHash(byteBuffer, frame, j, codeIncrement, length, firstCode);
              hashCallBidirectional(hash, reverseHash, s + startPos, internalId);
            } else {
              //System.err.println("hashCall seq=" + seq + " frame=" + frame + " j=" + j);
              hashCall(hash, internalId, s + startPos);
            }
          } //window
          end();
        }

        internalId++;
      }
      endSequence();
    }
    return totalLength;
  }

  private long reverseHash(final byte[] byteBuffer, final Frame frame, final int j, final int codeIncrement, final int length, final int firstCode) {
    mFunction.reset();
    long hash = 0L;
    for (int kk = 0, k = (j + mWindowSize - 1) * codeIncrement; kk < mWindowSize; kk++, k -= codeIncrement) {
      final byte b = frame.code(byteBuffer, length, k);
      final int c = b - firstCode;
      hash = mFunction.hashStep((byte) (3 - c));
    }
    return hash;
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

