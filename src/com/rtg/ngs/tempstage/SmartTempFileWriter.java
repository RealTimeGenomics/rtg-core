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
      throw new IllegalStateException("smart sam writer buffer distance (" + mBufferDistance + " bases, "
          + mRecordSet.size() + " records buffered) too small for record on ref " + record.getReferenceId()
          + " " + record.getStartPosition() + (record.isReverseStrand() ? "R" : "F")
          + " already written " + mLastWrittenRecord.getStartPosition() + (mLastWrittenRecord.isReverseStrand() ? "R" : "F")
          + " too small by " + (mLastWrittenRecord.getStartPosition() - record.getStartPosition() + 1));
    } else {
      Diagnostic.warning("SAM Record dropped due to excessive shift during alignment.");
    }
  }
}
