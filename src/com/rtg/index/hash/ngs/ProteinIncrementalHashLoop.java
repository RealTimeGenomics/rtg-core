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
package com.rtg.index.hash.ngs;

import java.io.IOException;

import com.rtg.index.hash.HashLoop;
import com.rtg.index.hash.ngs.protein.ProteinMask;
import com.rtg.launcher.ISequenceParams;
import com.rtg.launcher.HashingRegion;
import com.rtg.mode.Frame;
import com.rtg.mode.SequenceMode;
import com.rtg.mode.SequenceType;
import com.rtg.reader.SequencesReader;
import com.rtg.util.ProgramState;

/**
 * Loops over all hash codes derived from specified sequences in a reader.
 * This version increments the hash accumulator at each step.
 */
public abstract class ProteinIncrementalHashLoop extends HashLoop {

  protected final ProteinMask mProteinMask;

  /**
   * Constructs an instance and as a side effect makes calls to
   * <code>hashCall</code>.
   * @param windowSize number of codes to be included in a window used to make a hash.
   * @param stepSize number of steps to take per hash call
   * @param function used to construct the hash from the codes.
   */
  public ProteinIncrementalHashLoop(final int windowSize, final int stepSize, final ProteinMask function) {
    super(windowSize, stepSize);
    mProteinMask = function;
  }


  /**
   * Like other methods of the same name except without the extra information required.
   * (As it is provided in the constructor)
   * @param params parameters for the sequence
   * @return total number of codes read.
   * @throws IOException if an IO error occurs
   */
  public long execLoop(final ISequenceParams params) throws IOException {
    return execLoop(params, HashLoop.makeBuffer(params.reader()));
  }

  @Override
  public long execLoop(ISequenceParams params, final byte[] byteBuffer) throws IOException {
    return execLoop(byteBuffer, mWindowSize, mStepSize, params);
  }

  /**
   * @param byteBuffer buffer for data
   * @param windowSize should be <code>readLength - 1</code> (-1 to deal with missing last frame)
   * @param stepSize should be <code>readLength</code>
   * @param sequences the reader for data
   * @throws IOException If an I/O error occurs
   */
  private long execLoop(final byte[] byteBuffer, final int windowSize, final int stepSize, final ISequenceParams sequences) throws IOException {
    final SequencesReader reader = sequences.reader();
    final HashingRegion region = sequences.region();
    final long startSequence;
    final long endSequence;
    if (region != HashingRegion.NONE) {
      startSequence = region.getStart();
      endSequence = region.getEnd();
    } else {
      startSequence = 0;
      endSequence = reader.numberSequences();
    }
    final SequenceMode mode = sequences.mode();
    if (byteBuffer.length < reader.maxLength()) {
      throw new IllegalArgumentException("Allocated buffer too short. Allocated length=" + byteBuffer.length + " Required length=" + reader.maxLength());
    }
    //number of codes to be ignored (unknown etc).
    //final int firstCode = mode.codeType().firstValid();
    //the various possible frames (direction and phase for translation)
    final Frame[] frames = mode.allFrames();
    final int codeSize = mode.codeIncrement();
    final int posAdjustment = codeSize - windowSize * codeSize; //used to find the start of the read on the reference
    final int step = stepSize * codeSize;
    //Diagnostic.developerLog("stepSize: " + stepSize + " windowSize: " + windowSize + " codeSize: " + codeSize);
    //    final Timer readDelay = new Timer("Read_delay");
    //System.err.println("start=" + start + " end=" + end + " stepSize=" + stepSize + " step=" + step + " codeIncrement=" + codeIncrement);

    int internalId;
    if (reader.type() == SequenceType.DNA) {
      internalId = (region.getChunkId() * (int) reader.numberSequences() * frames.length) + (int) startSequence * frames.length;
    } else {
      internalId = (int) startSequence * frames.length;
    }
    //Diagnostic.developerLog("start: " + start + " internal:" + internalId);
    final int readLength = mProteinMask.readLength();
    //System.err.println(readLength);
    // To make sure we get all hits near the end of the template, we need to go the
    // read length less the window size extra positions during searching but not
    // during index building.  The SequenceType.DNA test below is a dubious way of
    // determining if we are in the build phase, since for map we build on DNA and
    // search on protein.
    final int overhang = reader.type() == SequenceType.DNA ? 0 : readLength - mProteinMask.windowSize();

    final long maxSequenceEver = reader.numberSequences();
    if (startSequence <= endSequence && startSequence < maxSequenceEver) {
      reader.seek(startSequence);
    }
    long totalLength = 0;
    for (long seq = startSequence; region.isInRange(seq) && seq < maxSequenceEver; seq++) {
      ProgramState.checkAbort();
      //System.err.println("seq=" + seq + " " +  reader.currentSequenceId() + " " + reader.getClass());
      final int startPos = region.getReferenceStart(seq, 0);
      final int endPos = region.getReferenceEnd(seq, 0, reader.currentLength());
      if (byteBuffer.length < endPos - startPos) { //silly programmer error
        throw new IllegalArgumentException("Allocated buffer too short. Allocated length=" + byteBuffer.length + " Required length=" + (endPos - startPos));
      }
      final int length = reader.readCurrent(byteBuffer, startPos, endPos - startPos);

      final int fullLength = reader.currentLength();
      totalLength += fullLength;
      //System.err.println("seq=" + seq + " " +  reader.currentSequenceId() + " " + reader.getClass());
      //      readDelay.start();


      templateSet(seq, length);
      nextSeq((int) (seq - startSequence), length);
      //      readDelay.stop(length);
      //System.out.println("chunk: " + region.getChunkId());
      for (final Frame frame : frames) {
        final int firstValid = frame.calculateFirstValid(startPos, length, fullLength);
        final int lastValid = frame.calculateLastValid(startPos, length, fullLength);
        //System.err.println("frame=" + frame);
        //System.out.println("frame=" + frame);
        next(seq, frame);
        final int phase = frame.phase();
        //final int limit = length - codeSize + 1 - phase + overhang;
        final int limit = lastValid - firstValid - codeSize + 1 - phase + overhang; // last valid protein

        //System.err.println("limit=" + limit + " phase=" + phase + " length=" + length + " codeIncrement=" + codeIncrement);
        mProteinMask.reset();
        for (int j = 0; j < limit; j += codeSize) { //iterate over entire reference (up to limit)
          //System.err.println("j=" + j);
          final byte b = frame.code(byteBuffer, length, j, firstValid, lastValid); //this internally deals with the code size for simplicity
          mProteinMask.hashStep(b);
          //System.out.print(Integer.toHexString(b & 0xff));
          final int pos = j + posAdjustment; //the 'actual start' of the read were this the last hash - used to deal with overhangs off the template.
          //System.err.println("seq=" + seq + " frame=" + frame + " j=" + j + " jcp=" + jcp + " step=" + step);
          //Diagnostic.developerLog("limit: " + limit + " isValid: " + mProteinMask.isValid() + " pos: " + pos + " j: " + j + " step: " + step + " pos % step: " + (pos % step));
          if (mProteinMask.isValid() && pos % step == 0) {
            assert pos >= mProteinMask.windowSize() - readLength - 1 : "jcp muther: " + pos + " windowSize: " + mProteinMask.windowSize() + " readLength: " + readLength + "   step: " + step + " jcpadj: " + posAdjustment + " codeinc: " + codeSize + " win: " + windowSize;  // -1 to account for the missing last frame
            //System.err.println("seq=" + seq + " frame=" + frame + " j=" + j + " step=" + step);

            hashCall(internalId, pos + readLength - 1); // add the readlength back to make this the last hash position, -1 due to the missing last frame
          }
        } //window
        //System.out.println();
        internalId++;
        end();
      }
      //System.err.println("Finished loop");
      reader.nextSequence();
      endSequence();

    }
    endAll();
    return totalLength;
  }

  @Override
  public void hashCall(final long hash, final int internalId, final int stepPosition) {
    throw new UnsupportedOperationException("Don't have a hash to supply here");
  }

  protected void templateSet(long name, int length) {
    // default do nothing
  }

  /**
   * As in {@link HashLoop#hashCall(long, int, int) } except that the hash is
   * stored by the internal hash function object
   * @param internalId the <code>internalId</code> of the sequence this hash is on
   * @param endPosition the end position of the window
   * @throws IOException If an I/O error occurs
   */
  public abstract void hashCall(final int internalId, final int endPosition) throws IOException;


  @Override
  public void next(final long seq, final Frame frame) {
    //default do nothing
  }

  @Override
  public void end() {
    //default do nothing
  }

  @Override
  public void endSequence() {
    //default do nothing
  }

  @Override
  public void endAll() throws IOException {
    //default do nothing
  }
}

