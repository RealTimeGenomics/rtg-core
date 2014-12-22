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

import java.io.Closeable;
import java.io.File;
import java.io.IOException;

import com.rtg.mode.SequenceMode;
import com.rtg.reader.SequencesReader;
import com.rtg.util.Params;
import com.rtg.util.intervals.LongRange;


/**
 */
public abstract class ReaderParams implements Params, Closeable {

  /**
   * Get the mode.
   * @return the mode.
   */
  public abstract SequenceMode mode();

  /**
   * Check if reader is closed.
   * @return true iff the reader is currently closed.
   */
  @Override
  public abstract boolean closed();

  /**
   * Get a SequencesReader for this sequence.
   * @return a SequencesReader for this sequence. A single reader is returned on succesive calls.
   */
  public abstract SequencesReader reader();

  /**
   * Get the lengths of the sequences in the reader.
   * @return the lengths.
   * @throws IOException if an I/O Error occurs
   */
  public abstract int[] lengths() throws IOException;

  /**
   * Get the length of the longest sequence.
   * @return the length of the longest sequence.
   */
  public abstract long maxLength();

  /**
   * If necessary carefully close the reader.
   * @throws IOException if an IO error occurs
   */
  @Override
  public abstract void close() throws IOException;

  /**
   * Returns a directory
   * @return directory
   */
  public abstract File directory();

  /**
   * Find if the region had to be adjusted to accommodate sloppy ending
   * @return if end of region sloppy adjusted, return new region, false otherwise
   */
  public abstract LongRange adjustedRegion();

}

