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
package com.rtg.util;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.Serializable;

import com.rtg.util.io.FileUtils;

/**
 * Handle creating streams in a directory.
 */
public class ResultStreamHandler implements Serializable {

  private final File mDir;
  /** If this is null, use <code>mDir</code> instead. */
  private final File mTempFilesDir;
  private final boolean mZip;

  /**
   * Constructor.
   * @param dir directory files are to be created in
   * @param tempFilesDir directory where intermediate files will be created
   * @param zip whether files should be compressed
   */
  public ResultStreamHandler(File dir, File tempFilesDir, boolean zip) {
    mDir = dir;
    mTempFilesDir = tempFilesDir;
    mZip = zip;
  }

  /**
   * Constructor.
   * @param dir directory files are to be created in
   * @param zip whether files should be compressed
   */
  public ResultStreamHandler(File dir, boolean zip) {
    this(dir, null, zip);
  }

  /**
   * Get the name of a child file in the output directory where all results are placed.
   * @param child the name of the child.
   * @return the name of the file.
   */
  public File file(final String child) {
    return new File(mDir, child);
  }

  /**
   * Gets the directory that will be used to store temporary files and
   * temporary subdirectories.
   *
   * @return a <code>File</code> representing the directory used to store temporary files.
   */
  public File tempDir() {
    return (mTempFilesDir != null) ? mTempFilesDir : mDir;
  }

  /**
   * Get the name of an intermediate file in the <code>tempdir</code> directory.
   * If no <code>tempdir</code> directory is set, the file will be put in the output directory.
   *
   * @param child the name of the child.
   * @return the name of the file.
   * @throws IOException if the file could not be determined
   */
  public File tempFile(final String child) throws IOException {
    final File dir = tempDir();
    if (!dir.exists() && !dir.mkdirs()) {
      throw new IOException("Unable to create temporary directory: " + dir.getPath());
    }
    return File.createTempFile(child + ".", "", dir);
  }

  /**
   * Create a file with the given name in the output directory, creating the directory
   * and any parent directory as needed. Also if output is compressed the compressed
   * suffix is appended to the file name.
   * @param name name of file
   * @return output stream into file
   * @throws IOException if there is a problem creating the stream
   */
  public OutputStream createFileStream(final String name) throws IOException {
    if (!mDir.exists() && !mDir.mkdirs()) {
      throw new IOException("Unable to create output file directory: " + mDir.getPath());
    }
    return FileUtils.createOutputStream(file(mZip ? name + FileUtils.GZ_SUFFIX : name));
  }
}
