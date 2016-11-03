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

package com.rtg.sam;

import java.io.IOException;

import com.rtg.util.ReorderingQueue;
import com.rtg.util.diagnostic.Diagnostic;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordCoordinateComparator;


/**
 * Class to reorder SAM records during output. Note that only one record per position is supported.
 */
public class SmartSamWriter extends ReorderingQueue<SAMRecord> {

  private static final int BUFFER_SIZE = 10000; // Assume records will be out of order by at most this amount.

  final SAMFileWriter mOut;

  /**
   * Constructor
   * @param out where records should eventually be sent
   */
  public SmartSamWriter(SAMFileWriter out) {
    super(BUFFER_SIZE, new SAMRecordCoordinateComparator() {
      public int compare(final SAMRecord samRecord1, final SAMRecord samRecord2) {
        final int cmp = super.compare(samRecord1, samRecord2);
        if (cmp != 0) {
          return cmp;
        }
        final SAMReadGroupRecord rg1 = samRecord1.getReadGroup();
        final SAMReadGroupRecord rg2 = samRecord2.getReadGroup();
        if (rg1 == rg2) {
          return 0;
        }
        return rg1 == null ? -1 : rg2 == null ? 1 : rg1.getReadGroupId().compareTo(rg2.getReadGroupId());
      }
    });
    mOut = out;
  }

  @Override
  protected String getReferenceName(SAMRecord record) {
    return record.getReferenceName();
  }

  @Override
  protected int getPosition(SAMRecord record) {
    return record.getAlignmentStart();
  }

  @Override
  protected void flushRecord(SAMRecord rec) throws IOException {
    mOut.addAlignment(rec);
  }

  @Override
  protected void reportReorderingFailure(SAMRecord rec) {
    Diagnostic.warning("SAMRecord dropped due to excessive out-of-order processing.\n" + rec.toString());
  }

  @Override
  public void close() throws IOException {
    super.close();
    mOut.close();
  }
}
