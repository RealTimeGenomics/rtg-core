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
package com.rtg.util.io;

import java.io.ByteArrayOutputStream;
import java.io.Closeable;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.PrintStream;

/**
 * For all those times you need to pass something to a print stream and then
 * dump the contents to a string
 */
public class MemoryPrintStream implements Closeable {

  private final ByteArrayOutputStream mBytes;
  private final PrintStream mPrint;
  private final LineWriter mWriter;

  /**
   * Construct an empty one
   */
  public MemoryPrintStream() {
    mBytes = new ByteArrayOutputStream();
    mPrint = new PrintStream(mBytes);
    mWriter = new LineWriter(new OutputStreamWriter(mBytes));
  }

  /**
   * Close the print stream.
   */
  @Override
  public void close() {
    mPrint.close();
  }

  /**
   * Retrieve the underlying ByteArrayOutputStream
   * @return the raw stream
   */
  public ByteArrayOutputStream outputStream() {
    return mBytes;
  }

  /**
   * Retrieve the PrintStream
   * @return the stream
   */
  public PrintStream printStream() {
    return mPrint;
  }

  /**
   * Retrieve the LineWriter
   * @return the stream
   */
  public LineWriter lineWriter() {
    return mWriter;
  }

  /**
   * Flushes first
   * @see java.io.ByteArrayOutputStream#toString()
   * @return string representation of contents of output
   */
  @Override
  public String toString() {
    mPrint.flush();
    flushWriter();
    return mBytes.toString();
  }

  private void flushWriter() {
    try {
      mWriter.flush();
    } catch (IOException e) {
      throw new RuntimeException(e); // Can't happen
    }
  }

  /**
   * @return contents of stream
   */
  public byte[] toByteArray() {
    mPrint.flush();
    flushWriter();
    return mBytes.toByteArray();
  }

  /**
   * Reset the underlying ByteArrayOutputStream
   * @see java.io.ByteArrayOutputStream#toString()
   */
  public void reset() {
    mBytes.reset();
  }
}
