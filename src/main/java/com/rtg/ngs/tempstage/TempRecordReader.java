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
