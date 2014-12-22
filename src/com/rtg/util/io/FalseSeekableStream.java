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

/**
 * Wraps a normal {@link InputStream} in a {@link SeekableStream} but does not
 * support seeking. Should not be passed around, instead only use inside a class
 * aware of its limitations.
 */
public class FalseSeekableStream extends SeekableStream {

  private final InputStream mStream;

  /**
   * Constructor
   * @param stream the stream to wrap
   */
  public FalseSeekableStream(InputStream stream) {
    mStream = stream;
  }

  /**
   * throws {@link UnsupportedOperationException}
   * @return nothing
   */
  @Override
  public long getPosition() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  /**
   * throws {@link UnsupportedOperationException}
   * @return nothing
   */
  @Override
  public long length() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  /**
   * throws {@link UnsupportedOperationException}
   * @param pos ignored
   */
  @Override
  public void seek(long pos) {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public int read() throws IOException {
    return mStream.read();
  }

  /**
   * throws {@link UnsupportedOperationException}
   * @return nothing
   */
  @Override
  public int available() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  /**
   */
  @Override
  public void close() throws IOException {
    mStream.close();
  }

  /**
   * throws {@link UnsupportedOperationException}
   * @param readlimit ignored
   */
  @Override
  public synchronized void mark(int readlimit) {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  /**
   * throws {@link UnsupportedOperationException}
   * @return nothing
   */
  @Override
  public boolean markSupported() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public int read(byte[] b) throws IOException {
    return mStream.read(b);
  }

  @Override
  public int read(byte[] b, int off, int len) throws IOException {
    return mStream.read(b, off, len);
  }

  /**
   * throws {@link UnsupportedOperationException}
   */
  @Override
  public synchronized void reset() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  /**
   * throws {@link UnsupportedOperationException}
   * @param n ignored
   * @return nothing
   */
  @Override
  public long skip(long n) {
    throw new UnsupportedOperationException("Not supported yet.");
  }
}
