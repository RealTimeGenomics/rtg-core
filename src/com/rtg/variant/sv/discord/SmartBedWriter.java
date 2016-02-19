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

package com.rtg.variant.sv.discord;

import java.io.IOException;
import java.io.OutputStream;
import java.io.Serializable;
import java.util.Comparator;

import com.rtg.util.ByteUtils;
import com.rtg.util.ReorderingQueue;
import com.rtg.util.diagnostic.Diagnostic;


/**
 * Class to reorder Bed records during output. Note that only one record per position is supported.
 */
public class SmartBedWriter extends ReorderingQueue<DiscordBedRecord> {

  private static final int BUFFER_SIZE = 10000; // Assume records will be out of order by at most this amount.

  /**
   * Comparator for out preferred bed output order
   */
  public static class BedPositionalComparator implements Comparator<DiscordBedRecord>, Serializable {
    @Override
    public int compare(DiscordBedRecord o1, DiscordBedRecord o2) {
      int res = o1.getSequenceName().compareTo(o2.getSequenceName());
      if (res == 0) {
        if (o1.getStart() > o2.getStart()) {
          res = 1;
        } else if (o1.getStart() < o2.getStart()) {
          res = -1;
        } else {
          res = Integer.compare(System.identityHashCode(o1), System.identityHashCode(o2)); // Ensure that none get clobbered
        }
      }
      return res;
    }
  }

  final OutputStream mOut;

  SmartBedWriter(OutputStream out) {
    super(BUFFER_SIZE, new BedPositionalComparator());
    mOut = out;
    //mOut.write(header.toString().getBytes()); // No header handling at present
  }

  @Override
  protected String getReferenceName(DiscordBedRecord record) {
    return record.getSequenceName();
  }

  @Override
  protected int getPosition(DiscordBedRecord record) {
    return record.getStart();
  }

  @Override
  protected void flushRecord(DiscordBedRecord rec) throws IOException {
    if (rec.isFiltered()) {
      mOut.write('#');
    }
    mOut.write(rec.toString().getBytes());
    ByteUtils.writeLn(mOut);
  }

  @Override
  protected void reportReorderingFailure(DiscordBedRecord rec) {
    Diagnostic.warning("DiscordBedRecord dropped due to excessive out-of-order processing.\n" + rec.toString());
  }

  @Override
  public void close() throws IOException {
    super.close();
    mOut.close();
  }
}
