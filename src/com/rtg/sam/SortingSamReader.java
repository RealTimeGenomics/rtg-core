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

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.Iterator;

import com.rtg.util.io.FileUtils;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

/**
 * Reads an input SAM/BAM in the given sort order by first writing out a temporary file
 */
public class SortingSamReader implements Iterable<SAMRecord>, Closeable {

  private final File mTmpFile;
  private final SamReader mReader;

  SortingSamReader(InputStream input, SAMFileHeader.SortOrder order, File tempDir) throws IOException {
    final SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault();
    final SamReader reader = samReaderFactory.open(SamInputResource.of(input));
    if (!reader.getFileHeader().getSortOrder().equals(order)) {
      try {
        final SAMFileHeader header = reader.getFileHeader();
        header.setSortOrder(order);
        mTmpFile = File.createTempFile("SortingSamReader", ".bam", tempDir);
        try (SAMFileWriter samFileWriter = new SAMFileWriterFactory().setTempDirectory(tempDir).makeBAMWriter(header, false, mTmpFile)) {
          for (SAMRecord rec : reader) {
            samFileWriter.addAlignment(rec);
          }
        }
      } finally {
        reader.close();
      }
      mReader = samReaderFactory.open(mTmpFile);
    } else {
      mTmpFile = null;
      mReader = reader;
    }

  }

  /**
   * Uses the system temporary directory to write the temporary file
   * @param samFile the SAM/BAM file to read
   * @param order the order to return the records in
   * @return the reader
   * @throws IOException if an IO error occurs
   */
  public static SortingSamReader reader(File samFile, SAMFileHeader.SortOrder order) throws IOException {
    return new SortingSamReader(FileUtils.createInputStream(samFile, true), order, null);
  }

  /**
   * @param samFile the SAM/BAM file to read
   * @param order the order to return the records in
   * @param tmpDir the directory in which to write the temporary file
   * @return the reader
   * @throws IOException if an IO error occurs
   */
  public static SortingSamReader reader(File samFile, SAMFileHeader.SortOrder order, File tmpDir) throws IOException {
    return new SortingSamReader(FileUtils.createInputStream(samFile, true), order, tmpDir);
  }

  /**
   * @return header to read SAM/BAM file
   */
  public SAMFileHeader getFileHeader() {
    return mReader.getFileHeader();
  }

  @Override
  public void close() throws IOException {
    IOException caught = null;
    boolean deleted;
    try {
      mReader.close();
    } catch (IOException exception) {
      caught = exception;
    } finally {
      deleted = mTmpFile.delete();
    }
    if (deleted) {
      if (caught != null) {
        throw caught;
      }
    } else {
      final IOException throwing = new IOException("Failed to delete temporary file: " + mTmpFile.getPath());
      if (caught == null) {
        throw throwing;
      } else {
        throwing.addSuppressed(caught);
        throw throwing;
      }
    }
  }

  @Override
  public Iterator<SAMRecord> iterator() {
    return mReader.iterator();
  }
}
