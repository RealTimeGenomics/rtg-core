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
import com.rtg.reference.Sex;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.Params;


/**
 */
public interface ISequenceParams extends Params {

  /**
   * Create a new version of this object whose start and end positions lie within the current range.
   * @param region the region to create sub params for
   * @return the new <code>SequenceParams</code>
   */
  ISequenceParams subSequence(HashingRegion region);

  /**
   * Get the mode.
   * @return the mode.
   */
  SequenceMode mode();

  /**
   * Get the directory.
   * @return the directory.
   */
  File directory();

  /**
   * Get the region that should be processed during query.
   * @return the region to be processed.
   */
  HashingRegion region();

  /**
   * What sex should this sequence be considered
   * @return the sex of the sequence
   */
  Sex sex();

  /**
   * Get the range of sequences that the reader should be restricted to
   * @return the region to be processed.
   */
  LongRange readerRestriction();

  /**
   * Get the total number of sequences.
   * @return the total number of sequences.
   */
  long numberSequences();

  /**
   * Get a SequencesReader for this sequence.
   * @return a SequencesReader for this sequence. A single reader is returned on successive calls.
   */
  SequencesReader reader();

  /**
   *  Get the reader parameters that specify the directory.
   * @return the reader parameters that specify the directory.
   */
  ReaderParams readerParams();

  /**
   * Get the length of the longest sequence.
   * @return the length of the longest sequence.
   */
  long maxLength();

  /**
   * If necessary close the reader.
   * @throws IOException If an IO error occurs
   */
  @Override
  void close() throws IOException;

}
