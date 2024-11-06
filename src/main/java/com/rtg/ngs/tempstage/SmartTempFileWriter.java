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


import java.io.IOException;
import java.io.OutputStream;
import java.util.Comparator;

import com.rtg.util.License;
import com.rtg.util.ReorderingQueue;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * Our smart temp file writer
 */
public class SmartTempFileWriter extends ReorderingQueue<BinaryTempFileRecord> {

  private final TempRecordWriter mOutputStream;

  /**
   * Makes our smart binary record temp file writer
   * @param out the data output stream to write to
   * @param comp comparator to sort records by
   * @param bufferDistance the initial buffer distance;
   */
  public SmartTempFileWriter(OutputStream out, Comparator<BinaryTempFileRecord> comp, int bufferDistance) {
    super(bufferDistance, comp);
    mOutputStream = new TempRecordWriterNio(out);
  }

  /**
   * Like the <code>addAlignment</code> method except will not write the record
   * if it already exists (assuming roughly template order of addition).
   * @param record record to add
   * @return true if record was added and will be written
   * @throws IOException if an IO exception occurs
   */
  public boolean addAlignmentHandleDuplicates(BinaryTempFileRecord record) throws IOException {
    return addRecord(record);
  }

  @Override
  public void close() throws IOException {
    super.close();
    mOutputStream.close();
  }

  @Override
  protected String getReferenceName(BinaryTempFileRecord record) {
    return Integer.toString(record.getReferenceId());
  }

  @Override
  protected int getPosition(BinaryTempFileRecord record) {
    return record.getStartPosition();
  }

  @Override
  protected void flushRecord(BinaryTempFileRecord rec) throws IOException {
    mOutputStream.writeRecord(rec);
  }

  @Override
  protected void reportReorderingFailure(BinaryTempFileRecord record) {
    // If we get here, it means either the aligner has shifted the
    // start position by more than it is allowed to, or the hits
    // from the index were more out of order than they should be!
    if (License.isDeveloper()) {
      for (final BinaryTempFileRecord brec : mRecordSet) {
        Diagnostic.developerLog("Rec: " + brec);
      }
      Diagnostic.developerLog("LastWrittenRecord: " + mLastWrittenRecord);
      Diagnostic.developerLog("Record out of order: " + record);
      throw new IllegalStateException("smart sam writer buffer distance (" + getWindowSize() + " bases, "
          + mRecordSet.size() + " records buffered) too small for record on ref " + record.getReferenceId()
          + " " + record.getStartPosition() + (record.isReverseStrand() ? "R" : "F")
          + " already written " + mLastWrittenRecord.getStartPosition() + (mLastWrittenRecord.isReverseStrand() ? "R" : "F")
          + " too small by " + (mLastWrittenRecord.getStartPosition() - record.getStartPosition() + 1));
    } else {
      Diagnostic.warning("SAM Record dropped due to excessive shift during alignment.");
    }
  }
}
