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
import java.io.RandomAccessFile;

/**
 * An extension of RandomAccessFile that reads and buffers large chucks of file content in order to reduce the number of disk reads/seeks.
 *
 */
public class BufferedRandomAccessFile extends RandomAccessFile {
  private static final int DEFAULT_BUFFER_SIZE = FileUtils.BUFFERED_STREAM_SIZE;

  private final byte[] mBuffer;
  private int mBufferEnd;
  private int mBufferPos;
  private long mRealPos;

  //private int mFillBufferCount = 0;
  //private int mReadBufCount = 0;
  //private int mReadCount = 0;
  //private int mSeekCount = 0;

  /**
   * Creates a random access file stream to read from, and optionally to write to, the file specified by the file argument.
   *
   * @param file the file
   * @param mode the access mode
   * @param bufsize size of internal buffer
   * @exception IOException if mode is readable and file does not exist.
   */
  public BufferedRandomAccessFile(final File file, final String mode, final int bufsize) throws IOException {
    super(file, mode);

    if (!"r".equals(mode)) {
      super.close();
      throw new IllegalArgumentException("Can only handle read only mode.");
    }
    if (bufsize < 1) {
      super.close();
      throw new IllegalArgumentException("Buffer size must be >= 1.");
    }
    mBuffer = new byte[bufsize];

    invalidate();
  }

  /**
   * Creates a random access file stream to read from, and optionally to write to, the file specified by the file argument.
   *
   * @param file the file
   * @param mode the access mode
   * @exception IOException if mode is readable and file does not exist.
   */
  public BufferedRandomAccessFile(final File file, final String mode) throws IOException {
    this(file, mode, DEFAULT_BUFFER_SIZE);
  }

  @Override
  public final int read() throws IOException {
    //mReadCount++;
    if (mBufferPos >= mBufferEnd) {
      if (fillBuffer() < 0) {
        return -1;
      }
    }
    return mBuffer[mBufferPos++] & 0xff;
  }

  private int fillBuffer() throws IOException {
    //mFillBufferCount++;
    final int n = super.read(mBuffer, 0, mBuffer.length);
    if (n >= 0) {
      mRealPos += n;
      mBufferEnd = n;
      mBufferPos = 0;
    }
    return n;
  }

  private void invalidate() throws IOException {
    mBufferEnd = 0;
    mBufferPos = 0;
    mRealPos = super.getFilePointer();
  }

  @Override
  public int read(final byte[] b, final int off, final int len) throws IOException {
    //mReadBufCount++;
    final int leftover = mBufferEnd - mBufferPos;
    if (len <= leftover) {
      System.arraycopy(mBuffer, mBufferPos, b, off, len);
      mBufferPos += len;
      return len;
    }
    //System.err.print(off + " " + len + " : " + leftover + " = " + mBufferEnd + " - " + mBufferPos);
    //System.err.println("\tRead");
    for (int i = 0; i < len; i++) {
      final int c = this.read();
      if (c != -1) {
        b[off + i] = (byte) c;
      } else if (i == 0) {
        return -1;
      } else {
        return i;
      }
    }
    return len;
  }

  @Override
  public long getFilePointer() {
    return mRealPos - mBufferEnd + mBufferPos;
  }

  @Override
  public void seek(final long pos) throws IOException {
    //mSeekCount++;
    final int n = (int) (mRealPos - pos);
    if (n >= 0 && n <= mBufferEnd) {
      mBufferPos = mBufferEnd - n;
    } else {
      super.seek(pos);
      invalidate();
    }
  }
  /*
  public void close() {
    super.close();
    System.err.println("SEEKS: " + mSeekCount);
    System.err.println("READS: " + mReadCount);
    System.err.println("READBUFS: " + mReadBufCount);
    System.err.println("FILLS: " + mFillBufferCount);

  }
  */
}

