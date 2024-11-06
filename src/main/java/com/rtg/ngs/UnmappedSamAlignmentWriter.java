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

package com.rtg.ngs;

import java.io.File;
import java.io.OutputStream;

import com.rtg.util.diagnostic.Diagnostic;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;

/**
 * The alignment writer for unmapped sam files
 */
public class UnmappedSamAlignmentWriter implements AutoCloseable {

  private static final int MAX_RECORDS_IN_RAM = 5000000;

  private final SAMFileHeader mHeader;
  private final SAMFileWriterFactory mFactory;
  private SAMFileWriter mSamUnmappedWriter = null;

  /**
   * Construct an unmapped record writer that can include unmated augmentation information
   *
   * @param tempFileDirectory directory to spill sorting data to disk in
   * @param header header to use for writing to SAM/BAM, if requested
   */
  public UnmappedSamAlignmentWriter(File tempFileDirectory, SAMFileHeader header) {
    mHeader = header;
    mFactory = new SAMFileWriterFactory();
    mFactory.setTempDirectory(tempFileDirectory);
    mFactory.setMaxRecordsInRam(MAX_RECORDS_IN_RAM);
  }

  /**
   * Tells the alignment writer where unmapped results should be written.
   * Must be called before any calls to <code>unmappedOut</code>
   *
   * @param unmappedOut output stream to send unmapped results to.
   * @param bam true if BAM should be written instead of SAM
   * @param presorted true if records do not require sorting
   * @param writeHeader true if output file should contain a header
   */
  public void initialiseUnmapped(OutputStream unmappedOut, boolean bam, boolean presorted, boolean writeHeader) {
    Diagnostic.userLog("Writing unmapped records");
    if (bam) {
      mSamUnmappedWriter = mFactory.makeBAMWriter(mHeader, presorted, unmappedOut, writeHeader, true, true);
    } else {
      mSamUnmappedWriter = mFactory.makeSAMWriter(mHeader, presorted, unmappedOut, writeHeader);
    }
  }

  /**
   * Write an unmapped record.
   * @param rec record to write
   */
  public void unmappedRecord(SAMRecord rec) {
    mSamUnmappedWriter.addAlignment(rec);
  }

  /**
   * Closes unmapped output
   */
  @Override
  public void close() {
    if (mSamUnmappedWriter != null) {
      mSamUnmappedWriter.close();
      mSamUnmappedWriter = null;
    }

  }
}
