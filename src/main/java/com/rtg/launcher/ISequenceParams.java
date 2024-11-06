/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
