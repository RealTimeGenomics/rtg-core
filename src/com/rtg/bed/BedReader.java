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

package com.rtg.bed;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;

import com.rtg.launcher.CommonFlags;
import com.rtg.tabix.BrLineReader;
import com.rtg.tabix.LineReader;
import com.rtg.tabix.TabixIndexer;
import com.rtg.tabix.TabixLineReader;
import com.rtg.util.intervals.RegionRestriction;
import com.rtg.util.io.FileUtils;

import net.sf.samtools.util.BlockCompressedInputStream;

/**
 * BED file reader class.
 */
public class BedReader implements Closeable {

  private final LineReader mIn;
  private final BedHeader mHeader;
  private BedRecord mCurrent;

  private final int mMinAnnotations;

  /**
   * create a new VCF reader.
   *
   * @param in where to read from
   * @throws IOException when IO or format errors occur.
   */
  public BedReader(BufferedReader in) throws IOException {
    this(in, 0);
  }

  /**
   * create a new BED reader with a minimum number of annotations
   *
   * @param in where to read from
   * @param minAnnotations the minimum number of annotations to require on each line of a file
   * @throws IOException when IO or format errors occur.
   */
  public BedReader(BufferedReader in, int minAnnotations) throws IOException {
    mMinAnnotations = minAnnotations;
    mIn = new BrLineReader(in);
    mHeader = parseHeader(mIn);
  }

  /**
   * Read BedRecords from a region of a block-compressed file
   * @param reader source of record lines
   * @param bedFile file reader was generated from
   * @throws IOException if an IO error occurs
   */
  public BedReader(TabixLineReader reader, File bedFile) throws IOException {
    mMinAnnotations = 0;
    mIn = reader;
    try (BrLineReader headerReader = new BrLineReader(new BufferedReader(new InputStreamReader(new BlockCompressedInputStream(bedFile))))) {
      mHeader = parseHeader(headerReader);
    }
    setNext();
  }

  /**
   * Open a <code>bed</code> reader, optionally with a region restriction
   * @param f <code>BED</code> file, optionally gzipped
   * @param region region to restrict to, null if whole file should be used
   * @return the reader
   * @throws IOException if an IO Error occurs
   */
  public static BedReader openBedReader(File f, RegionRestriction region) throws IOException {
    final boolean stdin = CommonFlags.isStdio(f);
    final BedReader bedr;
    if (region != null) {
      if (stdin) {
        throw new IOException("Cannot apply region restriction when reading from stdin");
      }
      bedr = new BedReader(new TabixLineReader(f, TabixIndexer.indexFileName(f), region), f);
    } else {
      bedr = new BedReader(new BufferedReader(new InputStreamReader(stdin ? System.in : FileUtils.createInputStream(f, true))));
    }
    return bedr;
  }

  private BedHeader parseHeader(LineReader in) throws IOException {
    final BedHeader header = new BedHeader();
    String line;
    while ((line = in.readLine()) != null) {
      try {
        if (line.startsWith("#") || line.startsWith("track") || line.startsWith("browser") || line.length() == 0) {
          header.addHeaderLine(line);
        } else {
          setCurrent(line);
          break;
        }
      } catch (final IllegalArgumentException e) {
        //illegal argument here means badly formed header
        throw new IOException(e.getMessage(), e);
      }
    }
    return header;
  }

  /**
   * Check if there is another record to get.
   * @return boolean true if there is another record to get
   */
  public boolean hasNext() {
    return mCurrent != null;
  }

  /**
   * Read the next record, if any.
   * @return the next record.
   * @throws IOException when IO or format errors occur.
   */
  public BedRecord next() throws IOException {
    if (mCurrent == null) {
      throw new IllegalStateException("No more records");
    }
    final BedRecord result = mCurrent;
    setNext();
    return result;
  }

  /**
   * Read the next record, if any.
   * @return true if there is a valid next record.
   * @throws IOException when IO or format errors occur.
   */
  private boolean setNext() throws IOException {
    String line;
    //TODO this would be problematic if there was ever a reference sequence named track
    while ((line = mIn.readLine()) != null && (line.length() == 0 || line.startsWith("track\t"))) { }
    if (line == null) {
      mCurrent = null;
      return false;
    }
    setCurrent(line);
    return true;
  }

  private void setCurrent(String line) throws IOException {
    try {
      mCurrent = BedRecord.fromString(line);
      if (mCurrent.getAnnotations().length < mMinAnnotations) {
        throw new IllegalArgumentException("Must have at least " + mMinAnnotations + " annotations");
      }
    } catch (final NumberFormatException e) {
      throw new IOException("Invalid BED line: Could not parse coordinates, line:" + line);
    } catch (final IllegalArgumentException | ArrayIndexOutOfBoundsException e) {
      //illegal argument == badly formed record
      throw new IOException("Invalid BED line: " + e.getMessage() + ", line:" + line, e);
    }
  }

  /**
   * Get the header.
   * @return the header.
   */
  public BedHeader getHeader() {
    return mHeader;
  }

  @Override
  public void close() throws IOException {
    mIn.close();
  }
}
