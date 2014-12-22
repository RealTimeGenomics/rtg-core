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
 * Stream that reads file contents into a buffer and then closes file handle until
 * it requires more data. File should only be open during read calls.
 */
public final class ClosedFileInputStream extends net.sf.samtools.util.SeekableStream {

  private static final int DEFAULT_BUFFER_SIZE = 1024 * 1024;

  private final File mFile;
  private final byte[] mBuffer;
  private final long mFileSize;
  private int mBufferPos;
  private int mBufferInUse;
  private long mFilePointer;
  private long mDiskSeeksDone = 0;

  /**
   * Constructor
   * @param file file to read from
   */
  public ClosedFileInputStream(File file) {
    this(file, DEFAULT_BUFFER_SIZE);
  }

  /**
   * @param file file to read from
   * @param bufferSize size of buffer
   */
  ClosedFileInputStream(File file, int bufferSize) {
    mBuffer = new byte[bufferSize];
    mFile = file;
    mFileSize = file.length();
  }

  /**
   * amount that can be read before opening the file again
   * @return the amount
   */
  @Override
  public int available() {
    return mBufferInUse - mBufferPos;
  }

  /**
   * does nothing since we don't keep the file open
   */
  @Override
  public void close() {
  }

  /**
   * @return true if we've reached end of file
   */
  @Override
  public boolean eof() {
    return mFilePointer == mFileSize;
  }

  /**
   * @return the files absolute path
   */
  @Override
  public String getSource() {
    return mFile.getAbsolutePath();
  }

  long getDiskSeeksDone() {
    return mDiskSeeksDone;
  }

  /**
   * @return current position within the file
   */
  public long position() {
    return mFilePointer;
  }

  @Override
  public long length() {
    return mFileSize;
  }

  @Override
  public void seek(long position) {
    final long bstart = mFilePointer - mBufferPos; // file position equiv to start of buffer
    final long bend = bstart + mBufferInUse;       // file position equiv to end of buffer

    if ((position >= bstart) && (position < bend)) {  // Allow free seek within the buffer without having to re-read
      mFilePointer = position;
      mBufferPos = (int) (position - bstart);
    } else {
      mFilePointer = position;
      mBufferPos = 0;
      mBufferInUse = 0;
    }
  }

  @Override
  public int read(byte[] buffer, int offset, int length) throws IOException {
    if (eof()) {
      return -1;
    }
    if (mBufferPos == mBufferInUse) {
      fillBuffer();
    }
    final int amount = Math.min(length, mBufferInUse - mBufferPos);
    System.arraycopy(mBuffer, mBufferPos, buffer, offset, amount);
    mBufferPos += amount;
    mFilePointer += amount;
    return amount;
  }

  @Override
  public int read() throws IOException {
    if (eof()) {
      return -1;
    }
    if (mBufferPos == mBufferInUse) {
      fillBuffer();
    }
    mFilePointer++;
    return mBuffer[mBufferPos++] & 0xff;
  }

  private void fillBuffer() throws IOException {
    mBufferPos = 0;
    mBufferInUse = 0;
    try (RandomAccessFile raf = new RandomAccessFile(mFile, "r")) {
      raf.seek(mFilePointer);
      mDiskSeeksDone++;
      int len;
      while (mBufferInUse != mBuffer.length && (len = raf.read(mBuffer, mBufferInUse, mBuffer.length - mBufferInUse)) != -1) {
        mBufferInUse += len;
      }
    }
  }
}
