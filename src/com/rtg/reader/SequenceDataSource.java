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
package com.rtg.reader;

import java.io.Closeable;
import java.io.IOException;

import com.rtg.mode.SequenceType;

/**
 * Provides access to a set of sequences.
 *
 */
public interface SequenceDataSource extends Closeable {

  /**
   * Get the type of the sequences (all sequences have the same type (DNA/Protein)).
   * @return the type of the sequences (non-null).
   */
  SequenceType type();

  /**
   * @return If reader contains quality data <code>true</code>, if not <code>false</code>
   */
  boolean hasQualityData();

  /**
   * Move to the next sequence.
   * @return true if there is a valid next sequence.
   * @throws IOException If in I/O error occurs
   */
  boolean nextSequence() throws IOException;

  /**
   * Get the length of the current sequence.
   * @return the length of the current sequence (&gt; 0).
   * @throws IllegalStateException if <code>nextSequence()</code> returned false on its last call.
   */
  int currentLength() throws IllegalStateException;

  /**
   * Get the name of the current sequence.
   * Will never be null and will always be 1 or more characters in length.
   * The set of characters that can occur in the name will be restricted to the
   * ASCII numbers 32 to 126 inclusive.
   * @return the name of the current sequence.
   * @throws IllegalStateException if <code>nextSequence()</code> returned false on its last call.
   * @throws IOException If in I/O error occurs
   */
  String name() throws IllegalStateException, IOException;

  /**
   * Returns the sequence data for the current sequence. This returns the internal byte array of the implementor. The array will only be filled up to <code>currentLength</code>.
   * @return the current sequence data.
   * @throws IllegalStateException if <code>nextSequence()</code> returned false on its last call.
   * @throws IOException If in I/O error occurs
   */
  byte[] sequenceData() throws IllegalStateException, IOException;

  /**
   * Returns the quality data for the current sequence. This returns the internal byte array of the implementor. The array will only be filled up to <code>currentLength</code>.
   * @return the current quality data.
   * @throws IllegalStateException if <code>nextSequence()</code> returned false on its last call.
   * @throws IOException If in I/O error occurs
   */
  byte[] qualityData() throws IllegalStateException, IOException;

  /**
   * Closes the data source.
   * @throws IOException If an I/O error occurs
   */
  @Override
  void close() throws IOException;

  /**
   * Enables dusting on the output sequence.
   * @param val True - enables dusting, False - disables it
   */
  void setDusting(boolean val);

  /**
   * Gets the number of dusted input residues.
   * @return the number of dusted input residues.
   */
  long getDusted();

  /**
   * Get the maximum sequence length.
   * @return the maximum sequence length.
   */
  long getMaxLength();

  /**
   * Get the minimum sequence length.
   * @return the minimum sequence length.
   */
  long getMinLength();

  /**
   * Return the total number of warnings that occurred
   * @return warning count
   */
  long getWarningCount();

}
