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

import java.io.IOException;

import com.rtg.sam.CigarFormatter;

/**
 * Class for presenting SAM/BAM format files to indexer
 */
class SamPositionReader extends AbstractPositionReader {

  private static final int NUM_COLUMNS = 10;
  private static final int REF_NAME_COLUMN = 2;
  private static final int START_POS_COLUMN = 3;
  private static final int CIGAR_COLUMN = 5;
  private static final int SEQ_COLUMN = 9;

  /**
   * Constructor
   * @param reader source of SAM file
   * @param skip number of lines to skip at start of file
   * @throws IOException if an IO error occurs.
   */
  SamPositionReader(BlockCompressedLineReader reader, int skip) throws IOException {
    super(reader, NUM_COLUMNS, '@', skip);
  }

  @Override
  protected void setReferenceName() throws IOException {
    mReferenceName = getColumn(REF_NAME_COLUMN);
  }

  @Override
  protected void setStartAndLength() throws IOException {
    mStartPosition = Integer.parseInt(getColumn(START_POS_COLUMN)) - 1;
    if (mStartPosition < 0) {
      mStartPosition = 0;
    }
    final String cigar = getColumn(CIGAR_COLUMN);
    if (cigar.equals("*")) {
      // For unmapped reads use length of sequence
      mLengthOnReference = getColumn(SEQ_COLUMN).length();
    } else {
      mLengthOnReference = CigarFormatter.cigarRefLength(cigar);
    }
  }


}
