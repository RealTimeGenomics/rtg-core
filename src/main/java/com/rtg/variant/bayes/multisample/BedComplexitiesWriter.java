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
package com.rtg.variant.bayes.multisample;

import java.io.IOException;
import java.io.OutputStream;

import com.rtg.util.ByteUtils;
import com.rtg.variant.bayes.multisample.ComplexRegion.RegionType;

/**
 * Write complexities for a single sequence into an output stream in BED format.
 */
public class BedComplexitiesWriter {

  private final OutputStream mOut;
  private final byte[] mSequenceName;
  private int mLastComplexitiesEnd;
  private int mStart = -1;
  private int mEnd = -1;
  private RegionType mType = null;

  /**
   * Construct an object for writing complex regions for a specific sequence to a BED file.
   * @param os output stream to write to
   * @param sequenceName sequence name
   * @param initialStartPosition position to start from
   */
  public BedComplexitiesWriter(OutputStream os, String sequenceName, int initialStartPosition) {
    mOut = os;
    mSequenceName = sequenceName.getBytes();
    mLastComplexitiesEnd = initialStartPosition;
  }

  private void writeRegion() throws IOException {
    if (mStart != -1) {
      mOut.write(mSequenceName);
      mOut.write(ByteUtils.TAB_BYTE);
      com.rtg.util.Utils.intWrite(mOut, mStart);
      mOut.write(ByteUtils.TAB_BYTE);
      com.rtg.util.Utils.intWrite(mOut, mEnd);
      mOut.write(ByteUtils.TAB_BYTE);
      mOut.write(mType.description());
      ByteUtils.writeLn(mOut);
    }
  }

  /**
   * Write complexities regions BED file.
   * @param complex complexities to write
   * @throws IOException if an IO error occurs.
   */
  public void write(final Complexities complex) throws IOException {
    if (complex.startOfChunk() != mLastComplexitiesEnd) {
      throw new IllegalStateException(complex.startOfChunk() + " != " + mLastComplexitiesEnd);
    }
    mLastComplexitiesEnd = complex.endOfChunk();
    for (final ComplexRegion cr : complex) {
      if (cr.getEnd() <= complex.startOfChunk() && cr.getStart() != cr.getEnd()) {
        continue;
      }
      if (cr.type() == mType && cr.getStart() <= mEnd + 1) { // + 1 to handle adjacent regions, end is exclusive
        // Combining regions of the same type.
        mEnd = Math.max(mEnd, cr.getEnd());
      } else {
        // cr is of a new type or is not immediately adjacent to saved region so write saved region and start a new one
        writeRegion();
        mStart = cr.getStart();
        mEnd = cr.getEnd();
        mType = cr.type();
      }
    }
  }

  /**
   * We are done writing complexities for current sequence.
   * @throws IOException if an IO error occurs.
   */
  public void finish() throws IOException {
    // Write any remaining region
    writeRegion();
  }
}
