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
import java.io.IOException;
import java.io.OutputStream;

import htsjdk.samtools.util.BlockCompressedOutputStream;

/**
 * An output stream that does the GZIP compression in a separate thread.
 * It also does buffering, before and after the GZIP compression and
 * within the pipe used to connect the GZIP thread to this thread.
 *
 * This class is final, because the helper thread starts in the
 * constructor.
 *
 */
public final class GzipAsynchOutputStream extends AsynchOutputStream {

  /** Whether to use bgzip instead of gzip */
  public static final boolean BGZIP = true; //Boolean.valueOf(System.getProperty("rtg.global-bgzip", "true"));

  /** Modern disks have about 64 Kb blocks, so make it that size? */
  public static final int DEFAULT_GZIP_BUFFER_SIZE = 64 * 1024;

  private static final int DEFAULT_GZIP_LEVEL = AdjustableGZIPOutputStream.DEFAULT_GZIP_LEVEL;

  // This is so we can do some parameter checking before calling super()
  private static OutputStream makeOutputStream(File file, int gzipSize) throws IOException  {
    if (file == null) {
      throw new IllegalArgumentException("File cannot be null");
    }
    if (!BGZIP) {
      return new BufferedOutputStreamFix(new AdjustableGZIPOutputStream(FileUtils.createOutputStream(file, false, false), gzipSize, DEFAULT_GZIP_LEVEL));
    } else {
    // Use BlockCompressedOutputStream so that the file is tabix compatible.
      return new BufferedOutputStreamFix(new BlockCompressedOutputStream(FileUtils.createOutputStream(file, false, false), null, DEFAULT_GZIP_LEVEL, false));
    }
  }

  /**
   * Create an asynchronous GZIP output stream with a pipe size of
   * <code>DEFAULT_BUFFER_SIZE</code> and a GZIP buffer size of
   * <code>DEFAULT_GZIP_BUFFER_SIZE</code>.
   *
   * @param file the output file.
   * @throws IOException on IO error.
   */
  public GzipAsynchOutputStream(File file) throws IOException {
    this(file, DEFAULT_PIPE_SIZE, DEFAULT_GZIP_BUFFER_SIZE);
  }

  /**
   * Create an asynchronous GZIP output stream to write to the given file.
   *
   * @param file the output file
   * @param pipeSize the size of the buffer between the threads.  At least 1 Kb.
   * @param gzipSize the buffer size of the decompression object.  At least 1 Kb.
   * @throws IOException on IO error.
   */
  GzipAsynchOutputStream(File file, int pipeSize, int gzipSize) throws IOException {
    super(makeOutputStream(file, gzipSize), pipeSize);
  }

 /**
  * Create an asynchronous GZIP output stream with a pipe size of
  * <code>DEFAULT_BUFFER_SIZE</code> and a GZIP buffer size of
  * <code>DEFAULT_GZIP_BUFFER_SIZE</code>.
  *
  * @param stream the output stream.
  * @param terminate true if the file should be terminated (for block compressed files)
  * @throws IOException on IO error.
  */
 public GzipAsynchOutputStream(OutputStream stream, boolean terminate) throws IOException {
   this(stream, DEFAULT_PIPE_SIZE, DEFAULT_GZIP_BUFFER_SIZE, terminate);
 }

 /**
  * Create an asynchronous GZIP output stream to write to the given stream.
  *
  * @param stream the output stream
  * @param pipeSize the size of the buffer between the threads.  At least 1 Kb.
  * @param gzipSize the buffer size of the decompression object.  At least 1 Kb.
  * @param terminated true to terminate the file (if block compressed)
  * @throws IOException on IO error.
  */
  public GzipAsynchOutputStream(OutputStream stream, int pipeSize, int gzipSize, boolean terminated) throws IOException {
    super(BGZIP ? new BlockCompressedOutputStream(stream, null, DEFAULT_GZIP_LEVEL, terminated) : new AdjustableGZIPOutputStream(stream), pipeSize);
    // Use BlockCompressedOutputStream so that the file is tabix compatible.
    //    super(new BlockCompressedOutputStream(stream, DEFAULT_GZIP_LEVEL), pipeSize, gzipSize);
  }
}

