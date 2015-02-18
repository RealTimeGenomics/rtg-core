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

/**
 * Position reader for <code>VCF</code> files, used for TABIX index support
 */
class VcfPositionReader extends AbstractPositionReader {

  private static final int NUM_COLUMNS = 4;
  private static final int REF_NAME_COLUMN = 0;
  private static final int START_POS_COLUMN = 1;
  private static final int REF_COLUMN = 3;

  /**
   * Constructor
   * @param reader source of SAM file
   * @param skip number of lines to skip at beginning of file
   * @throws IOException if an IO error occurs.
   */
  VcfPositionReader(BlockCompressedLineReader reader, int skip) throws IOException {
    super(reader, NUM_COLUMNS, '#', skip);
  }

  @Override
  protected void setReferenceName() throws IOException {
    mReferenceName = getColumn(REF_NAME_COLUMN);
  }

  @Override
  protected void setStartAndLength() throws IOException {
    mStartPosition = getIntColumn(START_POS_COLUMN) - 1;
    if (mStartPosition < 0) {
      mStartPosition = 0;
    }
    final String ref = getColumn(REF_COLUMN);

    mLengthOnReference = ref.length();
  }

}
