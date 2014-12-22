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
package com.rtg.util;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.Serializable;

import com.rtg.util.io.FileUtils;

/**
 * Handle creating streams in a directory
 */
public class ResultStreamHandler implements Serializable {

  private final File mDir;
  /** if this is null, use <code>mDir</code> instead. */
  private final File mTempFilesDir;
  private final boolean mZip;
  private final String mZipSuffix;

  /**
   * Constructor
   * @param dir Directory files are to be created in
   * @param tempFilesDir Directory where intermediate files will be created
   * @param zip Whether files should be compressed
   * @param zipSuffix suffix to use on compressed files
   */
  public ResultStreamHandler(File dir, File tempFilesDir, boolean zip, String zipSuffix) {
    mDir = dir;
    mTempFilesDir = tempFilesDir;
    mZip = zip;
    mZipSuffix = zipSuffix;
  }

  /**
   * Constructor
   * @param dir Directory files are to be created in
   * @param zip Whether files should be compressed
   * @param zipSuffix suffix to use on compressed files
   */
  public ResultStreamHandler(File dir, boolean zip, String zipSuffix) {
    this(dir, null, zip, zipSuffix);
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
    return FileUtils.createOutputStream(file(mZip ? name + mZipSuffix : name), mZip, false);
  }
}
