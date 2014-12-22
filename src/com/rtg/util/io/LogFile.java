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
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;

/**
 * A log backed by a file which can be deleted.
 */
public class LogFile implements LogStream {


  private final PrintStream mStream;

  private final File mFile;

  /**
   * @param file which is used to back the stream.
   */
  public LogFile(final File file) {
    if (file == null) {
      throw new NullPointerException();
    }
    mFile = file;
    FileOutputStream fileOutputStream;
    try {
      fileOutputStream = new FileOutputStream(file, true);
    } catch (final FileNotFoundException e) {
      fileOutputStream = null;
    }
    if (fileOutputStream == null) {
      mStream = null;
    } else {
      mStream = new PrintStream(fileOutputStream);
    }
  }

  @Override
  public void close() {
    if (mStream != null) {
      mStream.close();
    }
  }

  @Override
  public File file() {
    return mFile;
  }

  @Override
  public void removeLog() {
    close();
    if (mStream != null) {
      if (mFile.exists() && !mFile.delete()) {
        throw new RuntimeException("Unable to delete log file: " + mFile.getPath());
      }
    }
  }

  @Override
  public PrintStream stream() {
    return mStream;
  }

  @Override
  public String toString() {
    return "LogFile " + mFile.toString();
  }

}

