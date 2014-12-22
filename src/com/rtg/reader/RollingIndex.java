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
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.io.BufferedOutputStreamFix;
import com.rtg.util.io.FileUtils;

/**
 * A rolling index.
 *
 */
@TestClass(value = {"com.rtg.reader.SequencesWriterTest"})
class RollingIndex {
  DataOutputStream mOutput;
  private long mCurrent;
  private long mDataSize;
  private File mFile;

  public RollingIndex(final File file) {
    mCurrent = 0;
    mFile = file;
  }

  private void init() throws FileNotFoundException {
    mOutput = new DataOutputStream(new BufferedOutputStreamFix(new FileOutputStream(mFile), FileUtils.BUFFERED_STREAM_SIZE));
  }

  public void incrementCount() {
    mCurrent++;
  }

  public void incrementSize(long amount) {
    mDataSize += amount;
  }

  /**
   * Writes the count
   * @throws IOException if an I/O error occurs.
   */
  public void writeEntry() throws IOException {
    if (mOutput == null) {
      init();
    }
    mOutput.writeLong(mCurrent);
    mCurrent = 0;
    mOutput.writeLong(mDataSize);
    mDataSize = 0;
  }

  /**
   * Closes the output stream
   * @throws IOException if an I/O error occurs.
   */
  public void close() throws IOException {
    if (mOutput != null) {
      mOutput.close();
    }
  }
}
