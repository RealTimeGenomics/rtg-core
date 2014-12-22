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

import java.io.File;
import java.io.IOException;

import com.rtg.mode.SequenceType;

/**
 * Provides access to a set of sequences.
 *
 */
public interface SequencesReader extends AutoCloseable {
  /**
   * Get the type of the sequences (all sequences have the same type (DNA/Protein)).
   * @return the type of the sequences (non-null).
   */
  SequenceType type();

  /**
   * Get the sum of the lengths of all sequences.
   * @return the sum of the lengths of all sequences (&gt; 0).
   */
  long totalLength();

  /**
   * Get the length of the longest sequence.
   * @return the length of the longest sequence (&gt; 0).
   */
  long maxLength();

  /**
   * Get the length of the shortest sequence
   * @return the length of the shortest sequence. (&gt; 0)
   */
  long minLength();

  /**
   * Get the number of sequences.
   * @return the number of sequences (&gt; 0).
   */
  long numberSequences();

  /**
   * Returns the checksum for the sequence data
   * @return checksum
   */
  long dataChecksum();

  /**
   * Returns the checksum for the sequence qualities
   * @return checksum
   */
  long qualityChecksum();

  /**
   * Returns the checksum for the sequence names
   * @return checksum
   */
  long nameChecksum();

  /**
   * Returns the checksum for the sequence name suffixes
   * @return checksum
   */
  long suffixChecksum();

  /**
   * Return the residue counts for current preread
   * residues are stored on the basis of ordinal
   * @return residue counts array,
   */
  long[] residueCounts();

  /**
   * Return histogram of N's
   * @return array containing counts
   */
  long[] histogram();

  /**
   * Return position histogram of N's
   * @return array containing counts
   */
  long[] posHistogram();

  /**
   * Return quality scores average histogram
   * @return array containing counts
   */
  double globalQualityAverage();

  /**
   * Return the average quality per position
   * @return array containing counts
   */
  double[] positionQualityAverage();
  /**
   * Number of blocks of N's
   * @return the number
   */
  long nBlockCount();

  /**
   * Longest single N block
   * @return the length
   */
  long longestNBlock();

  /**
   * returns true if result from <code>histogram()</code> is valid
   * @return true if histogram available
   */
  boolean hasHistogram();

  /**
   * @return left and right
   */
  PrereadArm getArm();

  /**
   * @return sequencing technology
   */
  PrereadType getPrereadType();

  /**
   * Get a unique id which is the same on left and right arms of Complete Genomics
   * data. For use with older data only.
   * @return GUID
   */
  SdfId getSdfId();


  /**
   * Returns the SDF version of the underlying source if appropriate
   * @return the version of the SDF store, or -1 if it is not an SDF store
   */
  long sdfVersion();

  /**
   * Return the the path to any source file/directory relating to this reader.  If
   * this reader is not backed by files then it is permissible to return
   * null.
   *
   * @return directory
   */
  File path();

  /**
   * @return If reader contains quality data <code>true</code>, if not <code>false</code>
   */
  boolean hasQualityData();

  /**
   * @return If reader contains sequence names <code>true</code>, if not <code>false</code>
   */
  boolean hasNames();

  /// Sequence fetching methods

  /**
   * Position the reader at the specified sequence.
   * If you intend to call {@link SequencesReader#nextSequence()} after this method,
   * the <code>sequenceId</code> should be (<code>sequenceId - 1</code>).
   * @param sequenceId the sequence to seek to.
   * @throws IOException If in I/O error occurs
   */
  void seek(long sequenceId) throws IOException;

  /**
   * Move to the next sequence.
   * @return true if there is a valid next sequence.
   * @throws IOException If in I/O error occurs
   */
  boolean nextSequence() throws IOException;

  /**
   * Get the identifier for the current sequence.
   * Will initially be set to 0 and is incremented by <code>nextSequence()</code>.
   * @return the identifier for the current sequence ( &gt;= 0).
   * @throws IllegalStateException if <code>nextSequence()</code> returned false on its last call.
   */
  long currentSequenceId() throws IllegalStateException;

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
  String currentName() throws IllegalStateException, IOException;

  /**
   * Get the full name of the current sequence. Formally this is the result of the current suffix appended to the current name
   * Will never be null and will always be 1 or more characters in length.
   * @return the name of the current sequence.
   * @throws IllegalStateException if <code>nextSequence()</code> returned false on its last call.
   * @throws IOException If in I/O error occurs
   */
  String currentFullName() throws IllegalStateException, IOException;

  /**
   * Get the suffix of the current name, normally this will be everything after and including the first whitespace character
   * from the source data file.
   * @return the suffix
   * @throws IllegalStateException if <code>nextSequence()</code> returned false on its last call.
   * @throws IOException If in I/O error occurs
   */
  String currentNameSuffix() throws IllegalStateException, IOException;

  /**
   * Reads current sequence  into the supplied array.
   * @param dataOut array to read data into
   * @return length of sequence
   * @throws IllegalArgumentException If <code>dataOut</code> does not have enough length to store sequence.
   * @throws IllegalStateException if <code>nextSequence()</code> returned false on its last call.
   * @throws IOException If in I/O error occurs
   */
  int readCurrent(byte[] dataOut) throws IllegalArgumentException, IllegalStateException, IOException;

  /**
   * Reads sequence data into the supplied array.
   * @param dataOut array to read data into
   * @param start the start offset within the sequence to read from
   * @param length the number of residues to read
   * @return length of sequence read
   * @throws IllegalArgumentException If <code>dataOut</code> does not have enough length to store sequence.
   * @throws IOException If in I/O error occurs
   */
  int readCurrent(byte[] dataOut, int start, int length) throws IllegalArgumentException, IOException;

  /**
   * Reads current quality into the supplied array.
   * @param dest array to read data into
   * @return length of quality, 0 if <code>hasQualityData()</code> is false
   * @throws IllegalArgumentException If <code>dataOut</code> does not have enough length to store quality.
   * @throws IllegalStateException if <code>nextSequence()</code> returned false on its last call.
   * @throws IOException If in I/O error occurs
   */
  int readCurrentQuality(byte[] dest) throws IllegalArgumentException, IllegalStateException, IOException;

  /**
   * Reads current quality into the supplied array.
   * @param dest array to read data into
   * @param start the start offset within the sequence to read from
   * @param length the number of residues to read
   * @return length of quality, 0 if <code>hasQualityData()</code> is false
   * @throws IllegalArgumentException If <code>dataOut</code> does not have enough length to store quality.
   * @throws IllegalStateException if <code>nextSequence()</code> returned false on its last call.
   * @throws IOException If in I/O error occurs
   */
  int readCurrentQuality(byte[] dest, int start, int length) throws IllegalArgumentException, IllegalStateException, IOException;

  /**
   * Returns the length of the requested sequence
   * @param sequenceIndex the sequence id
   * @return the length of the requested sequence
   */
  int length(long sequenceIndex);

  /**
   * Get the name of the specified sequence.
   * Will never be null and will always be 1 or more characters in length.
   * The set of characters that can occur in the name will be restricted to the
   * ASCII numbers 32 to 126 inclusive.
   * @param sequenceIndex Sequence to read.
   * @return the name of the current sequence.
   * @throws IOException If in I/O error occurs
   */
  String name(long sequenceIndex) throws IOException;

  /**
   * Get the name suffix of the specified sequence.
   * @param sequenceIndex Sequence to read.
   * @return the name of the current sequence.
   * @throws IOException If in I/O error occurs
   */
  String nameSuffix(long sequenceIndex) throws IOException;

  /**
   * Get the full name of the specified sequence.
   * Will never be null and will always be 1 or more characters in length.
   * The set of characters that can occur in the name will be restricted to the
   * ASCII numbers 32 to 126 inclusive.
   * @param sequenceIndex Sequence to read.
   * @return the name of the current sequence.
   * @throws IOException If in I/O error occurs
   */
  String fullName(long sequenceIndex) throws IOException;

  /**
   * Reads sequence data into the supplied array.
   * @param sequenceIndex Sequence to read.
   * @param dataOut array to read data into
   * @return length of sequence
   * @throws IllegalArgumentException If <code>dataOut</code> does not have enough length to store sequence.
   * @throws IOException If in I/O error occurs
   */
  int read(long sequenceIndex, byte[] dataOut) throws IllegalArgumentException, IOException;

  /**
   * Reads sequence data into the supplied array.
   * @param sequenceIndex Sequence to read.
   * @param dataOut array to read data into
   * @param start the start offset within the sequence to read from
   * @param length the number of residues to read
   * @return length of sequence read
   * @throws IllegalArgumentException If <code>dataOut</code> does not have enough length to store sequence.
   * @throws IOException If in I/O error occurs
   */
  int read(long sequenceIndex, byte[] dataOut, int start, int length) throws IllegalArgumentException, IOException;

  /**
   * Reads quality data into the supplied array.
   * @param sequenceIndex Sequence to read quality for.
   * @param dest array to read data into
   * @return length of quality, 0 if <code>hasQualityData()</code> is false
   * @throws IllegalArgumentException If <code>dataOut</code> does not have enough length to store quality.
   * @throws IOException If in I/O error occurs
   */
  int readQuality(long sequenceIndex, byte[] dest) throws IllegalArgumentException, IOException;

  /**
   * Reads quality data into the supplied array.
   * @param sequenceIndex Sequence to read quality for.
   * @param dest array to read data into
   * @param start the start offset within the sequence to read from
   * @param length the number of quality values to read
   * @return length of quality, 0 if <code>hasQualityData()</code> is false
   * @throws IllegalArgumentException If <code>dataOut</code> does not have enough length to store quality.
   * @throws IOException If in I/O error occurs
   */
  int readQuality(long sequenceIndex, byte[] dest, int start, int length) throws IllegalArgumentException, IOException;

  /**
   * If appropriate ensures any backing file is closed.
   *
   * @throws IOException If in I/O error occurs
   */
  @Override
  void close() throws IOException;

  /**
   * Return an object which can be used to get names for sequences.
   * @return names
   * @throws IOException If in I/O error occurs
   */
  PrereadNamesInterface names() throws IOException;

  /**
   * count the number of residue in the sequences between <code>start</code>(inclusive) and <code>end</code>(exclusive).
   * @param start sequence id of first sequence.
   * @param end  sequence id of last sequence + 1
   * @return the length
   * @throws IOException If in I/O error occurs
   */
  long lengthBetween(long start, long end) throws IOException;

  /**
   * Return all the sequence lengths
   *
   * @param start starting sequence id
   * @param end ending sequence id (excl)
   * @return array of sequence lengths
   * @throws IOException If in I/O error occurs
   */
  int[] sequenceLengths(long start, long end) throws IOException;

  /**
   * @return another copy of this reader which can be used independently
   */
  SequencesReader copy();

  /**
   * Whether input was compressed, generally only required for filtering purposes
   * @return true if compressed
   */
  boolean compressed();

  /**
   * Reset the Sequences Reader to initial state.
   */
   void reset();

   /**
   * Get the contents of the read-me file as a string if this reader has a directory.
   * @return the contents of the read-me, or null if it does not exist.
   * @throws IOException if there is an error reading the file.
    */
   String getReadMe() throws IOException;
}
