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

package com.rtg.ngs;

import java.io.File;
import java.io.OutputStream;

import com.rtg.util.diagnostic.Diagnostic;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;

/**
 * The alignment writer for unmapped sam files
 */
public class UnmappedSamAlignmentWriter {

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
  public void close() {
    if (mSamUnmappedWriter != null) {
      mSamUnmappedWriter.close();
      mSamUnmappedWriter = null;
    }

  }
}
