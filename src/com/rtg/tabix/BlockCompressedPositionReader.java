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
package com.rtg.tabix;

import java.io.Closeable;
import java.io.IOException;
import java.util.List;

/**
 * Interface for presenting data to <code>BAM/TABIX</code> style indexers.
 */
public interface BlockCompressedPositionReader extends Closeable {

  /**
   * @return string representation of current record
   */
  String getRecord();

  /**
   * Seek to position defined by virtual offset
   * @param virtualOffset virtual offset as used by <code>TABIX/BAM</code> indexes
   * @throws IOException if an IO error occurs
   */
  void seek(long virtualOffset) throws IOException;

  /**
   * @return true if more records exist
   */
  boolean hasNext();

  /**
   * advances to the next record.
   * @throws IOException if an IO Error occurs
   */
  void next() throws IOException;

  /**
   * @return 0-based start position on the reference of the current record
   */
  int getStartPosition();

  /**
   * @return length on the reference of the current record
   */
  int getLengthOnReference();

  /**
   * @return name of reference sequence for current record
   */
  String getReferenceName();

  /**
   * @return an integer for the current reference name, should match the indices produced by {@link BlockCompressedPositionReader#getSequenceNames()}
   */
  int getReferenceId();

  /**
   * Note: virtual file offset defined as <code>offset_into_file_of_compressed_block &lt;&lt; 16 | offset_into_uncompressed_block</code>
   * @return the virtual file offset of the start of the current record.
   */
  long getVirtualOffset();

  /**
   * Note: virtual file offset defined as <code>offset_into_file_of_compressed_block &lt;&lt; 16 | offset_into_uncompressed_block</code>
   * @return the virtual file offset of the start of the next record.
   */
  long getNextVirtualOffset();

  /**
   * @return the bin number of the current record. See {@link TabixIndexer#reg2bin(int, int)}
   */
  int getBinNum();

  /**
   * @return return true if the current record is unmapped. Always false for <code>TABIX</code> style indexing.
   */
  boolean isUnmapped();

  /**
   * @return return true if the current record has a reference sequence. Always true for <code>TABIX</code> style indexing.
   */
  boolean hasReference();

  /**
   * @return return true if the current record has coordinates. Always true for <code>TABIX</code> style indexing.
   */
  boolean hasCoordinates();

  /**
   * @return List of sequence names in order they appear in the file being indexed.
   */
  List<String> getSequenceNames();
}
