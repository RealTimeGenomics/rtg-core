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

import java.io.Closeable;
import java.io.File;
import java.io.PrintStream;

/**
 * A wrapper for a stream to be used for logging.
 * Provides facilities for deletion and closing.
 * These may be ignored if the stream is not backed by a file.
 */
public interface LogStream extends Closeable {

  /**
   * Get the stream.
   * @return the stream.
   */
  PrintStream stream();

  /**
   * Close the underlying stream if this makes sense.
   * Otherwise the call does nothing.
   */
  @Override
  void close();

  /**
   * If the stream is backed by an underlying file then delete it.
   * Otherwise do nothing.
   */
  void removeLog();

  /**
   * Return a backing file if one exists.
   * @return the backing file (null if there isn't one).
   */
  File file();
}

