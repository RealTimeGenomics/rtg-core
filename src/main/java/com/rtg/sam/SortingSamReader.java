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
    if (reader.getFileHeader().getSortOrder() != order) {
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
    final boolean deleted;
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
      if (caught != null) {
        throwing.addSuppressed(caught);
      }
      throw throwing;
    }
  }

  @Override
  public Iterator<SAMRecord> iterator() {
    return mReader.iterator();
  }
}
