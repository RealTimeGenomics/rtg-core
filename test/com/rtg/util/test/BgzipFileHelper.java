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

package com.rtg.util.test;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

import com.rtg.util.Resources;
import com.rtg.util.io.FileUtils;

import net.sf.samtools.util.BlockCompressedOutputStream;

/**
 */
public final class BgzipFileHelper {

  private BgzipFileHelper() { }

  /**
   * Convenience method for turning a resource into a block compressed GZIP file that can be used as input.
   * @param resource name (classpath) of the resource.
   * @param file a <code>File</code> to write to
   * @return a <code>File</code> containing the string content (same as <code>file</code>).
   * @throws IOException if an error occurs
   */
  public static File resourceToBgzipFile(final String resource, final File file) throws IOException {
    try (InputStream stream = Resources.getResourceAsStream(resource)) {
      return streamToBgzipFile(stream, file);
    }
  }
  /**
   * Writes the contents of the given input <code>stream</code> to the given block compressed
   * gzipped file.
   *
   * @param stream an <code>InputStream</code>
   * @param file a <code>File</code> to write to
   * @return a <code>File</code> containing the contents of the stream
   * @exception IOException if an error occurs.
   * @exception NullPointerException if the stream is null
   */
  public static File streamToBgzipFile(final InputStream stream, final File file) throws IOException {
    if (stream == null) {
      throw new NullPointerException("null stream given");
    }
    if (file == null) {
      throw new NullPointerException("null file given");
    }
    try (OutputStream out = new BlockCompressedOutputStream(file)) {
      final byte[] b = new byte[FileUtils.BUFFER_SIZE];
      int len = stream.read(b);
      while (len > 0) {
        out.write(b, 0, len);
        len = stream.read(b);
      }
    }
    return file;
  }
  /**
   * Create a BGZIP file
   * @param data data to write
   * @param f file to write to
   * @return file written to
   * @throws IOException if an IO error occurs
   */
  public static File bytesToBgzipFile(byte[] data, File f) throws IOException {
    try (BlockCompressedOutputStream out = new BlockCompressedOutputStream(f)) {
      out.write(data);
    }
    return f;
  }
}
