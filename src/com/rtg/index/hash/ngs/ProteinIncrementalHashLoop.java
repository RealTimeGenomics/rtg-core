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
  protected final boolean mTemplate;

  /**
   * Constructs an instance and as a side effect makes calls to
   * <code>hashCall</code>.
   * @param windowSize number of codes to be included in a window used to make a hash.
   * @param stepSize number of steps to take per hash call
   * @param function used to construct the hash from the codes.
   * @param isTemplate true if the hashloop is used for the template, false when used for reads
   */
  public ProteinIncrementalHashLoop(final int windowSize, final int stepSize, final ProteinMask function, boolean isTemplate) {
    super(windowSize, stepSize);
    mProteinMask = function;
    mTemplate = isTemplate;
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
    return execLoop(params, byteBuffer, mWindowSize, mStepSize);
  }

  /**
   * @param sequences the reader for data
   * @param byteBuffer buffer for data
   * @param windowSize should be <code>readLength - 1</code> (-1 to deal with missing last frame)
   * @param stepSize should be <code>readLength</code>
   * @throws IOException If an I/O error occurs
   */
  private long execLoop(final ISequenceParams sequences, final byte[] byteBuffer, final int windowSize, final int stepSize) throws IOException {
    final SequencesReader reader = sequences.reader();
    final HashingRegion region = sequences.region();
    final long startSequence;
    if (region != HashingRegion.NONE) {
      startSequence = region.getStart();
    } else {
      startSequence = 0;
    }
    if (byteBuffer.length < reader.maxLength()) {
      throw new IllegalArgumentException("Allocated buffer too short. Allocated length=" + byteBuffer.length + " Required length=" + reader.maxLength());
    }

    final SequenceMode mode = sequences.mode();
    final Frame[] frames = mode.allFrames();
    final int codeSize = mode.codeIncrement();
    final int posAdjustment = codeSize - windowSize * codeSize; //used to find the start of the read on the reference
    final int step = stepSize * codeSize;

    int internalId;
    if (reader.type() == SequenceType.DNA) {
      internalId = (region.getChunkId() * (int) reader.numberSequences() * frames.length) + (int) startSequence * frames.length;
    } else {
      internalId = (int) startSequence * frames.length;
    }
    final int readLength = mProteinMask.readLength();
    // To make sure we get all hits near the end of the template, we need to go the
    // read length less the window size extra positions during template searching but not
    // during index building.
    final int overhang = mTemplate ? readLength - mProteinMask.windowSize() : 0;

    final long maxSequenceEver = reader.numberSequences();
    long totalLength = 0;
    for (long seq = startSequence; region.isInRange(seq) && seq < maxSequenceEver; ++seq) {
      ProgramState.checkAbort();
      final int fullLength = reader.length(seq);
      final int startPos = region.getReferenceStart(seq, 0);
      final int endPos = region.getReferenceEnd(seq, 0, fullLength);
      if (byteBuffer.length < endPos - startPos) { //silly programmer error
        throw new IllegalArgumentException("Allocated buffer too short. Allocated length=" + byteBuffer.length + " Required length=" + (endPos - startPos));
      }
      final int length = reader.read(seq, byteBuffer, startPos, endPos - startPos);
      totalLength += fullLength;

      templateSet(seq, length);
      nextSeq((int) (seq - startSequence), length);
      for (final Frame frame : frames) {
        final int firstValid = frame.calculateFirstValid(startPos, length, fullLength);
        final int lastValid = frame.calculateLastValid(startPos, length, fullLength);
        //System.out.println("frame=" + frame);
        next(seq, frame);
        final int phase = frame.phase();
        //final int limit = length - codeSize + 1 - phase + overhang;
        final int limit = lastValid - firstValid - codeSize + 1 - phase + overhang; // last valid protein

        mProteinMask.reset();
        for (int j = 0; j < limit; j += codeSize) { //iterate over entire reference (up to limit)
          final byte b = frame.code(byteBuffer, length, j, firstValid, lastValid); //this internally deals with the code size for simplicity
          mProteinMask.hashStep(b);
          final int pos = j + posAdjustment; //the 'actual start' of the read were this the last hash - used to deal with overhangs off the template.
          if (mProteinMask.isValid() && pos % step == 0) {
            assert pos >= mProteinMask.windowSize() - readLength - 1 : "jcp muther: " + pos + " windowSize: " + mProteinMask.windowSize() + " readLength: " + readLength + "   step: " + step + " jcpadj: " + posAdjustment + " codeinc: " + codeSize + " win: " + windowSize;  // -1 to account for the missing last frame
            hashCall(internalId, pos + readLength - 1); // add the readlength back to make this the last hash position, -1 due to the missing last frame
          }
        } //window
        ++internalId;
        end();
      }
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
  public void endAll() {
    //default do nothing
  }
}

