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
        res = Integer.compare(o1.getStart(), o2.getStart());
        if (res == 0) {
          res = Integer.compare(o1.getEnd(), o2.getEnd());
          if (res == 0) {
            res = Integer.compare(System.identityHashCode(o1), System.identityHashCode(o2)); // Ensure that none get clobbered
          }
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
    Diagnostic.warning("DiscordBedRecord dropped due to excessive out-of-order processing.\n" + rec);
  }

  @Override
  public void close() throws IOException {
    super.close();
    mOut.close();
  }
}
