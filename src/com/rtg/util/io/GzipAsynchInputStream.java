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
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;

import com.rtg.util.gzip.GzipUtils;


/**
 * An input stream that does the GZIP deflation in a separate thread.
 * It also does buffering, before and after the Gzip deflation and
 * within the pipe used to connect the GZIP thread to this thread.
 *
 * This class is final, because the helper thread starts in the
 * constructor.
 *
 */
public final class GzipAsynchInputStream extends AsynchInputStream {

  /** Modern disks have about 64 Kb blocks, so make it that size? */
  public static final int DEFAULT_GZIP_BUFFER_SIZE = 64 * 1024;

  // This is so we can do some parameter checking before calling super()
  private static InputStream makeInputStream(File file, int gzipSize) throws IOException {
    if (file == null) {
      throw new IllegalArgumentException("File cannot be null");
    }
    final FileInputStream is = new FileInputStream(file);
    try {
      return GzipUtils.createGzipInputStream(is, gzipSize);
    } catch (IOException e) {
      is.close();
      throw e;
    }
  }

  /**
   * Create an asynchronous GZIP input stream to read the given file.
   *
   * @param file the input file
   * @param pipeSize the size of the buffer between the threads.  At least 1 Kb.
   * @param gzipSize the buffer size of the decompression object.
   * @throws IOException on IO error.
   */
  public GzipAsynchInputStream(File file, int pipeSize, int gzipSize) throws IOException {
    super(makeInputStream(file, gzipSize), pipeSize, gzipSize);
  }

  /**
   * Create an asynchronous GZIP input stream with <code>DEFAULT_GZIP_BUFFER_SIZE</code>.
   *
   * @param file the input file.
   * @throws IOException on IO error.
   */
  public GzipAsynchInputStream(File file) throws IOException {
    this(file, DEFAULT_PIPE_SIZE, DEFAULT_GZIP_BUFFER_SIZE);
  }

  /**
   * Create an asynchronous GZIP input stream with <code>DEFAULT_GZIP_BUFFER_SIZE</code>.
   * @param is the input stream.
   * @throws IOException on IO error.
   */
  public GzipAsynchInputStream(InputStream is) throws IOException {
    super(GzipUtils.createGzipInputStream(is));
  }

}

