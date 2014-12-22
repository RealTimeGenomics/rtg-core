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
package com.rtg.reader;

import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;

import com.rtg.util.io.BufferedOutputStreamFix;
import com.rtg.util.io.FileUtils;


/**
 * Tracks pairs of data and pointers
 */
class NameFilePair {
  private final OutputStream mNameData;
  private final DataOutputStream mPointers;
  private final long mLimit;
  private int mDataSize;

  /**
   * Creates an instance
   * @param names Destination file for names
   * @param pointers Destination file for pointers
   * @param limit Max file size in bytes (currently ignored)
   * @throws IOException When IO errors occur
   */
  public NameFilePair(final File names, final File pointers, final long limit) throws IOException {
    mNameData = new BufferedOutputStreamFix(new FileOutputStream(names), FileUtils.BUFFERED_STREAM_SIZE);
    mPointers = new DataOutputStream(new BufferedOutputStreamFix(new FileOutputStream(pointers), FileUtils.BUFFERED_STREAM_SIZE));
    mLimit = limit;
  }

  /**
   * Checks whether writing name will break file size limit
   * @param length length of the name
   * @return true
   */
  public boolean canWriteName(final int length) {
    return mDataSize + length + 1 <= mLimit
      && mPointers.size() + 4 <= mLimit;
  }

  /**
   * Writes a name to files.
   * @param name Name to write
   * @throws IOException When an IO error occurs
   */
  public void writeName(final String name) throws IOException {
    mPointers.writeInt(mDataSize);
    final byte[] rawChars = name.getBytes();
    mNameData.write(rawChars);
    mNameData.write(0); //it worked for c (kinda)
    mDataSize += rawChars.length + 1;
  }

  public String forceWriteName(final String name) throws IOException {
    final long remaining = mLimit - mDataSize;
    if (remaining > name.length() + 1) {
      writeName(name);
      return name;
    } else {
      final String truncated = name.substring(0, (int) (remaining - 1));
      writeName(truncated);
      return truncated;
    }
  }

  /**
   * Closes internal OutputStreams
   * @throws IOException
   */
  public void close() throws IOException {
    mNameData.close();
    mPointers.close();
  }

  /**
   * @return number of ASCII values written including null characters
   */
  public long valuesWritten() {
    return mDataSize;
  }
}
