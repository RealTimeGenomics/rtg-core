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


/**
 * Stream manager for reading labels
 * should only be used by DefaultSequencesReader
 */
class LabelStreamManager extends AbstractStreamManager {
  /**
   * Creates a the stream manager
   * @param dir Directory containing sequence data
   * @param numberSequences Number of sequences in directory
   * @throws IOException If an I/O Error occurs
   */
  LabelStreamManager(final File dir, final long numberSequences, final String indexFile, final String dataPrefix, final String pointerPrefix, long dataIndexVersion, DataFileOpenerFactory openerFactory) throws IOException {
    super(dir, numberSequences, indexFile, dataPrefix, pointerPrefix, dataIndexVersion, openerFactory.getLabelOpener());
  }

  /**
   * Creates a stream manager for names
   */
  static LabelStreamManager getNameStreamManager(File dir, long numberSequences, long dataIndexVersion, DataFileOpenerFactory openerFactory) throws IOException {
    return new LabelStreamManager(dir, numberSequences, SdfFileUtils.LABEL_INDEX_FILENAME, SdfFileUtils.LABEL_DATA_FILENAME, SdfFileUtils.LABEL_POINTER_FILENAME, dataIndexVersion, openerFactory);
  }

  /**
   * Creates a stream manager for name suffixes
   */
  static LabelStreamManager getSuffixStreamManager(File dir, long numberSequences, long dataIndexVersion, DataFileOpenerFactory openerFactory) throws IOException {
    return new LabelStreamManager(dir, numberSequences, SdfFileUtils.LABEL_SUFFIX_INDEX_FILENAME, SdfFileUtils.LABEL_SUFFIX_DATA_FILENAME, SdfFileUtils.LABEL_SUFFIX_POINTER_FILENAME, dataIndexVersion, openerFactory);
  }

  @Override
  protected void seekImpl(final long seqNum) throws IOException {
    final long pointerpos = seqNum - mCurrentLower;
    mPointers.randomAccessFile().seek(pointerpos * 4L);

    final long seqpos = readInt(mPointers.randomAccessFile()); //.readInt();
    final long length;
    if (mPointers.randomAccessFile().length() - mPointers.randomAccessFile().getPosition() >= 4L) {
      int nextseqpos;
      nextseqpos = readInt(mPointers.randomAccessFile()); //.readInt();
      length = nextseqpos - seqpos;
    } else {
      length = mIndex.dataSize(mIndexedSequenceFileNumber) - seqpos; //mData.randomAccessFile().length() - seqpos;
    }
    if (length - 1 > Integer.MAX_VALUE || length < 0) {
      //we cast to int in various places
      throw new CorruptSdfException();
    }
    mDataLength = length - 1;
    mData.randomAccessFile().seek(seqpos);
  }
}

