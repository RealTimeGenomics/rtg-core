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

package com.rtg.util.store;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.Reader;

import com.rtg.util.integrity.IntegralAbstract;
import com.rtg.util.io.FileUtils;

/**
 */
public class StoreFileProxy extends IntegralAbstract implements StoreFile {

  private final File mFile;

  private final String mName;
  /**
   * Construct <code>StoreDirectory</code> rooted at the specified file.
   * @param file name of a file to be constructed.
   */
  StoreFileProxy(final File file, final String name) {
    mFile = file;
    mName = name;
  }

  @Override
  public OutputStream outputStream() throws FileNotFoundException {
    final FileOutputStream fos;
    fos = new FileOutputStream(mFile);
    return fos;
  }

  @Override
  public InputStream inputStream() throws IOException {
    final FileInputStream fos;
    fos = new FileInputStream(mFile);
    return fos;
  }

  @Override
  public String content() throws IOException {
    try (Reader reader = new InputStreamReader(inputStream())) {
      return FileUtils.readerToString(reader);
    }
  }

  @Override
  public String name() {
    return mName;
  }

  @Override
  public boolean integrity() {
    return true;
  }

}
