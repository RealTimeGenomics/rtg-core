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
import java.io.InputStream;

import com.reeltwo.jumble.annotations.TestClass;

/**
 * Class to get around the bizarre lack of reasonable interfaces on <code>java.io.RandomAccessFile</code>
 */
@TestClass("com.rtg.util.io.RandomAccessFileStreamTest")
public abstract class SeekableStream extends InputStream {

  /**
   * Sets the position offset, measured from the beginning of this
   * stream, at which the next read or write occurs.  The offset may be
   * set beyond the end of the stream. Setting the offset beyond the end
   * of the file does not change the underlying structure.  The length will
   * change only by writing after the offset has been set beyond the end
   * of the stream.
   *
   * @param      pos   the offset position, measured in bytes from the
   *                   beginning of the stream, at which to set the position.
   * @exception  IOException  if <code>pos</code> is less than
   *                          <code>0</code> or if an I/O error occurs.
   */
  public abstract void seek(long pos) throws IOException;

  /**
   * Returns the current offset in this stream.
   *
   * @return     the offset from the beginning of the stream, in bytes,
   *             at which the next read or write occurs.
   * @exception  IOException  if an I/O error occurs.
   */
  public abstract long getPosition() throws IOException;

  /**
   * Returns the length of this stream.
   *
   * @return     the length of this stream, measured in bytes.
   * @exception  IOException  if an I/O error occurs.
   */
  public abstract long length() throws IOException;

  /**
   */
  @Override
  public void close() throws IOException {
    super.close();
  }
}
