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

import java.io.Closeable;
import java.io.IOException;
import java.io.OutputStream;

import com.rtg.util.ByteUtils;


/**
 * Writer to write BED records out into a <code>.bed</code> output stream.
 *
 */
public class BedWriter implements Closeable {

  /** The character that indicates that the line is a comment */
  public static final byte COMMENT_CHAR = '#';

  private final OutputStream mOut;

  /**
   * create a new BED writer
   *
   * @param out stream to write to
   */
  public BedWriter(OutputStream out) {
    if (out == null) {
      throw new NullPointerException("output stream cannot be null");
    }
    mOut = out;
  }

  /**
   * Write out a BED record
   *
   * @param record record to write
   * @throws java.io.IOException if error
   */
  public void write(BedRecord record) throws IOException {
    writeLine(record.toString());
  }

  /**
   * Write a line as a comment (will automatically be prepended with comment char)
   * @throws java.io.IOException if error
   */
  protected void writeComment(String line) throws IOException {
    mOut.write(COMMENT_CHAR);
    writeLine(line);
  }

  /**
   * Write an extra line (e.g. track line) to output stream
   * @throws java.io.IOException if error
   */
  protected void writeLine(String line) throws IOException {
    mOut.write(line.getBytes());
    ByteUtils.writeLn(mOut);
  }

  @Override
  public void close() throws IOException {
    mOut.close();
  }

}
