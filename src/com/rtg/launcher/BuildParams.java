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
package com.rtg.launcher;

import java.io.File;
import java.io.IOException;

import com.rtg.mode.SequenceMode;
import com.rtg.reader.SequencesReader;
import com.rtg.util.Params;
import com.rtg.util.Utils;

/**
 * Holds a set of parameters that are sufficient to enable the initial creation
 * and build of an index.
 *
 */
public final class BuildParams implements Params {

  static long size(final int windowSize, final int stepSize, final ISequenceParams sequenceParams) throws IOException {
    final SequenceMode mode = sequenceParams.mode();
    final HashingRegion region = sequenceParams.region(); // XXX This should probably be sequenceParams.readerRestriction()??
    final long st;
    final long en;
    if (region != HashingRegion.NONE) {
      st = region.getStart();
      en = region.getEnd();
    } else {
      st = 0;
      en = sequenceParams.reader().numberSequences();
    }
    final SequencesReader reader = sequenceParams.reader();
    final long length = reader.lengthBetween(st, en);
    //Diagnostic.developerLog("start=" + st + " end=" + en + " size=" + length + " windowSize=" + windowSize + " stepSize=" + stepSize);
    return size(length, en - st, windowSize, stepSize, mode.numberFrames(), mode.codeIncrement());
  }

  /**
   * Compute the total number of hash windows given a sequence length and other values.
   * @param length of sequence.
   * @param numSequences number of sequences
   * @param windowSize window size.
   * @param stepSize step size.
   * @param frames number of frames.
   * @param codeIncrement number of residues in one code (needed for translation).
   * @return the total number of hash windows.
   */
  public static long size(final long length, long numSequences, final int windowSize, final int stepSize, final int frames, final int codeIncrement) {
    final long l = (length / codeIncrement) - windowSize;
    if (l < 0) {
      return 0;
    }
    long size = ((l / stepSize) + 1) * frames;
    if (stepSize > windowSize) { // When stepsize is larger than windowsize the above formula may underestimate the number of hashes per sequence by at most one hash per sequence
      size += numSequences;
    }
    //Diagnostic.developerLog("size=" + size + " length=" + length + " l=" + l + " windowSize=" + windowSize + " stepSize=" + stepSize + " codeIncrement=" + codeIncrement + " frames=" + frames);
    //System.err.println("size=" + size + " length=" + length + " l=" + l + " windowSize=" + windowSize + " stepSize=" + stepSize + " codeIncrement=" + codeIncrement);
    return size;
  }

  /**
   * calculate the total number of hash windows based on sequences window and step sizes
   * @return the number of hash windows
   * @throws IOException if an IO error occurs
   */
  public long calculateSize() throws IOException {
    return size(mWindowSize, mStepSize, mSequenceParams);
  }

  /**
   * Creates a BuildParams builder.
   * @return the builder.
   */
  public static BuildParamsBuilder builder() {
    return new BuildParamsBuilder();
  }


  /**
   * A builder class for <code>BuildParams</code>.
   */
  public static class BuildParamsBuilder {

    protected ISequenceParams mSequenceParams;

    protected int mWindowSize;

    protected int mStepSize;

    protected int mTypeBits; // Calculated automatically

    /**
     * create the builder
     */
    public BuildParamsBuilder() {
    }

    /**
     * Sets the sequences that will be used to build this index.
     * @param sequences the sequence parameters.
     * @return this builder, so calls can be chained.
     */
    public BuildParamsBuilder sequences(final ISequenceParams sequences) {
      mSequenceParams = sequences;
      return this;
    }

    /**
     * Sets window size.
     * @param windowSize number of residues in a hash window.
     * @return this builder, so calls can be chained.
     */
    public BuildParamsBuilder windowSize(final int windowSize) {
      mWindowSize = windowSize;
      return this;
    }

    /**
     * Sets number of residues between the start of succesive windows.
     * @param stepSize number of residues between the start of succesive windows.
     * @return this builder, so calls can be chained.
     */
    public BuildParamsBuilder stepSize(final int stepSize) {
      mStepSize = stepSize;
      return this;
    }

    /**
     * Creates a BuildParams using the current builder
     * configuration.
     * @return the new BuildParams
     */
    public BuildParams create() {
      if (mSequenceParams == null) {
        throw new NullPointerException("BuildParams requires SequenceParams to be set.");
      }
      mTypeBits = mSequenceParams.mode().codeType().bits();
      return new BuildParams(this);
    }
  }



  private final ISequenceParams mSequenceParams;

  /** Number of bits needed to represent a single residue (3 or 5). */
  private final int mTypeBits;

  private final int mWindowSize;

  private final int mStepSize;

  /**
   * Create a set of parameters to use in building a SLIM index.
   *
   * @param builder the builder object.
   */
  public BuildParams(final BuildParamsBuilder builder) {
    mSequenceParams = builder.mSequenceParams;
    mTypeBits = builder.mTypeBits;
    mWindowSize = builder.mWindowSize;
    mStepSize = builder.mStepSize;

    if (mSequenceParams == null) {
      throw new IllegalArgumentException();
    }
    if (mTypeBits != 2 && mTypeBits != 5) {
      throw new IllegalArgumentException();
    }
    if (mWindowSize < 1) {
      throw new IllegalArgumentException();
    }
    if (mStepSize < 1) {
      throw new IllegalArgumentException();
    }
  }

  /**
   * Create a builder with all the values set to those of this object.
   * @return a builder
   */
  public BuildParamsBuilder cloneBuilder() {
    return new BuildParamsBuilder().windowSize(windowSize()).stepSize(stepSize()).sequences(sequences());
  }

  /**
   * Create a new parameters object where everything is the same except for the range in the subsequence.
   * @param region of subsequence to construct sub params for
   * @return the new <code>BuildParams</code>.
   * @throws IOException If in I/O error occurs
   */
  public BuildParams subSequence(final HashingRegion region) throws IOException {
    //System.err.println(sequences.getClass().getName() + "=sequences class " + sequences + "=sequences");
    return cloneBuilder().sequences(sequences().subSequence(region)).create();
  }

  /**
   * Get the directory where the subject sequences are to be found.
   * @return the directory where the subject sequences are to be found.
   */
  public File directory() {
    return mSequenceParams.directory();
  }

  /**
   * Get the sequence parameters.
   * @return the sequence parameters.
   */
  public ISequenceParams sequences() {
    return mSequenceParams;
  }

  /**
   * Get the number of residues in a hash window.
   * @return the number of residues in a hash window.
   */
  public int windowSize() {
    return mWindowSize;
  }

  /**
   * Get the number of residues between the start of succesive windows.
   * @return the number of residues between the start of succesive windows.
   */
  public int stepSize() {
    return mStepSize;
  }

  @Override
  public void close() throws IOException {
    mSequenceParams.close();
  }

  @Override
  public int hashCode() {
    return Utils.pairHash(Utils.pairHash(mSequenceParams.hashCode(), mTypeBits), Utils.pairHash(mWindowSize, mStepSize));
  }

  @Override
  public boolean equals(final Object obj) {
    if (this == obj) {
      return true;
    }
    if (!(obj instanceof BuildParams)) {
      return false;
    }
    final BuildParams that = (BuildParams) obj;
    return mTypeBits == that.mTypeBits
      && mWindowSize == that.mWindowSize
      && mStepSize == that.mStepSize
      && mSequenceParams.equals(that.mSequenceParams);
  }

  @Override
  public String toString() {
    return " seq={" + mSequenceParams + "} "
     + " window=" + mWindowSize
     + " step=" + mStepSize;
  }
}
