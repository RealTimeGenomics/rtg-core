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

import java.io.File;
import java.io.PrintStream;

/**
 */
public class LogSimple implements LogStream {

  private final PrintStream mStream;

  /**
   * @param stream which is being wrapped by this class.
   */
  public LogSimple(final PrintStream stream) {
    if (stream == null) {
      throw new NullPointerException();
    }
    mStream = stream;
  }

  @Override
  public void close() {
    if (mStream != System.err) {
      mStream.close();
    }
  }

  @Override
  public void removeLog() {
    close();
  }

  @Override
  public PrintStream stream() {
    return mStream;
  }

  @Override
  public File file() {
    return null;
  }

  @Override
  public String toString() {
    return "LogSimple";
  }

}

