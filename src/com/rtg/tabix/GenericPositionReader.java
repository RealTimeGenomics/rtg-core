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
 * Reader for providing data for generic files to be <code>TABIX</code> indexed
 */
public class GenericPositionReader extends AbstractPositionReader {
  private final boolean mOneBased;
  private final int mStartCol;
  private final int mSeqCol;
  private final int mEndCol;

  /**
   * @param reader read positions out of this reader
   * @param options the options to use while reading positions
   * @throws IOException when IO breaks
   */
  public GenericPositionReader(BlockCompressedLineReader reader, TabixIndexer.TabixOptions options) throws IOException {
    super(reader, noColsNeeded(options), (char) options.mMeta, options.mSkip);
    mOneBased = !options.mZeroBased;
    mStartCol = options.mStartCol;
    mEndCol = options.mEndCol;
    mSeqCol = options.mSeqCol;
  }

  private static int noColsNeeded(TabixIndexer.TabixOptions options) {
    //cols are zero based in options class
    return Math.max(Math.max(options.mStartCol, options.mEndCol), options.mSeqCol) + 1;
  }

  @Override
  protected void setReferenceName() throws IOException {
    mReferenceName = getColumn(mSeqCol);
  }

  @Override
  protected void setStartAndLength() throws IOException {
    mStartPosition = Integer.parseInt(getColumn(mStartCol)) - (mOneBased ? 1 : 0);
    if (mEndCol != mStartCol) {
      mLengthOnReference = Integer.parseInt(getColumn(mEndCol)) - mStartPosition;
    } else {
      mLengthOnReference = 1;
    }
  }

}
