/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
package com.rtg.index.params;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

import com.rtg.launcher.OutputDirParams;
import com.rtg.util.Constants;
import com.rtg.util.ObjectParams;
import com.rtg.util.Utils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.Integrity;
import com.rtg.util.io.FileUtils;

/**
 * Holds enough information to generate a hit counts writer.
 */
public class CountParams extends ObjectParams implements OutputDirParams, Integrity {

  private final File mCountDir;

  private final int mTopN;

  private final int mMinHits;

  private final long mMaxFileSize;

  private final boolean mZip;

  /**
   * @param countDir the directory where the hit counts will be written.
   * @param topN the number of best hits to be used.
   * @param minHits ignore  hit unless it occurs this number of times.
   * @param maxFileSize the maximum file size allowed for individual files.
   * @param zip true iff outputs are to be zipped.
   */
  public CountParams(final File countDir, final int topN, final int minHits, final long maxFileSize, boolean zip) {
    mZip = zip;
    mTopN = topN;
    mMinHits = minHits;
    mCountDir = countDir;
    mMaxFileSize = maxFileSize;
    integrity();
  }

  /**
   * Create with a default maximum file size.
   * @param countDir the directory where the hit counts will be written.
   * @param topN the number of best hits to be used.
   * @param minHits ignore  hit unless it occurs this number of times.
   * @param zip true iff outputs are to be zipped.
   */
  public CountParams(final File countDir, final int topN, final int minHits, boolean zip) {
    this(countDir, topN, minHits, Constants.MAX_FILE_SIZE, zip);
  }

  /**
   * Get a stream to the output file.
   * @param name file name
   * @return the stream.
   * @throws IOException if an I/O error occurs.
   */
  public OutputStream outStream(final String name) throws IOException {
    if (!directory().isDirectory() && !directory().mkdirs()) {
      throw new IOException("Unable to create directory \"" + directory().getPath() + "\"");
    }
    final File filename = file(mZip ? name + FileUtils.GZ_SUFFIX : name);
    return FileUtils.createOutputStream(filename);
  }

  /**
   * Get the name of a child file in the output directory where all results are placed.
   * @param child the name of the child.
   * @return the name of the file.
   */
  public File file(final String child) {
    return new File(mCountDir, child);
  }

  /**
   * Get the directory.
   * @return the directory.
   */
  @Override
  public File directory() {
    return mCountDir;
  }

  /**
   * Get  the maximum number of hits for a single query sequence.
   * @return the maximum number of hits for a single query sequence.
   */
  public int topN() {
    return mTopN;
  }

  /**
   * Get  the minimum number of hits required for a hit to be reported.
   * @return the minimum number of hits required for a hit to be reported.
   */
  public int minHits() {
    return mMinHits;
  }

  /**
   * Get the maximum file size allowed for individual files.
   * @return the maximum file size allowed for individual files.
   */
  public long maxFileSize() {
    return mMaxFileSize;
  }
  @Override
  public int hashCode() {
    return Utils.pairHash(Utils.pairHash(mCountDir.hashCode(), (int) mMaxFileSize), Utils.pairHash(mTopN, mMinHits));
  }

  @Override
  public boolean equals(final Object obj) {
    if (this == obj) {
      return true;
    }
    if (obj == null) {
      return false;
    }
    final CountParams that = (CountParams) obj;
    return this.mCountDir.equals(that.mCountDir)
    && this.mTopN == that.mTopN
    && this.mMinHits == that.mMinHits
    && this.mMaxFileSize == that.mMaxFileSize;
  }

  @Override
  public String toString() {
    return " CountParams directory=" + mCountDir + " topN=" + mTopN + " minHits=" + mMinHits + " max. file size=" + mMaxFileSize;
  }

  @Override
  public final boolean integrity() {
    Exam.assertTrue(mCountDir != null);
    Exam.assertTrue(mTopN >= 1);
    Exam.assertTrue(mMinHits >= 1);
    Exam.assertTrue(mMaxFileSize >= Constants.MINIMUM_FILE_CHUNK_SIZE);
    return true;
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    return true;
  }

}

