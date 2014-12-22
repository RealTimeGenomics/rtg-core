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

import java.io.IOException;
import java.io.OutputStream;
import java.util.zip.Deflater;
import java.util.zip.GZIPOutputStream;

/**
 * An override of the standard <code>GZIPOutputStream</code> so that
 * we can get access to the deflation compression level.
 *
 *
 */
public class AdjustableGZIPOutputStream extends GZIPOutputStream {

  /** Modern disks have about 64 Kb blocks, so make it that size? */
  public static final int DEFAULT_GZIP_BUFFER_SIZE = 64 * 1024;

  /** Benchmarking shows that 2 gives best CPU to compression tradeoff */
  public static final int DEFAULT_GZIP_LEVEL = 2;


  /**
   * Creates a GZIP output stream with a given compression level.
   * The fastest level seems to be 2 or 3.
   * The compression gets significantly better at 6, but is almost
   * twice as slow as the fastest.
   * The default compression for the standard GZIPOutputStream is -1,
   * which seems to correspond to about level 6.
   *
   * @param out the stream to send the compressed data to.
   * @param gzipSize the size of the GZIP buffers.
   * @param level -1 or 0 to 9 where 0 means no compression and 9 is maximum.
   * @throws IOException on IO error.
   */
  public AdjustableGZIPOutputStream(OutputStream out, int gzipSize, int level) throws IOException {
    super(out, gzipSize);
    //Diagnostic.developerLog("new AdjustableGZIPOutputStream(" + gzipSize + ", " + level + ")");
    assert -1 <= level && level <= 9;
    def = new Deflater(level, true);
  }

  /**
   * Creates a GZIP output stream with default compression level and buffer size.
   *
   * @param out the stream to send the compressed data to.
   * @throws IOException on IO error.
   */
  public AdjustableGZIPOutputStream(OutputStream out) throws IOException {
    this(out, DEFAULT_GZIP_BUFFER_SIZE, DEFAULT_GZIP_LEVEL);
  }
}
