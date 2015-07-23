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
package com.rtg.ngs.tempstage;

import java.io.Closeable;
import java.io.IOException;

/**
 */
public interface TempRecordReader extends Closeable {

  /**
   * Reads a record from the input stream and inserts the data into a temp file record
   * @return the temp file record which was read, or null when the end of the stream is reached
   * @throws IOException if something horrible happens when reading a record
   */
  BinaryTempFileRecord readRecord() throws IOException;

  /**
   * Return a temp file record factory of the given type
   */
  class RecordFactory {
    private final boolean mPaired;
    private final boolean mLegacy;
    private final boolean mCg;
    private final boolean mUnfiltered;

    /**
     * @param paired are records paired end at sequencing
     * @param legacyCigars are we using legacy cigars
     * @param cg cg mode
     * @param unfiltered unfiltered mode
     */
    public RecordFactory(boolean paired, boolean legacyCigars, boolean cg, boolean unfiltered) {
      mPaired = paired;
      mLegacy = legacyCigars;
      mCg = cg;
      mUnfiltered = unfiltered;
    }

    BinaryTempFileRecord createRecord() {
      return new BinaryTempFileRecord(mPaired, mLegacy, mCg, mUnfiltered);
    }
  }
}
