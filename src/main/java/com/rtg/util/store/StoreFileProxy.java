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
   * @param name the name to associate with the file
   */
  StoreFileProxy(final File file, final String name) {
    mFile = file;
    mName = name;
  }

  @Override
  public OutputStream outputStream() throws FileNotFoundException {
    return new FileOutputStream(mFile);
  }

  @Override
  public InputStream inputStream() throws IOException {
    return new FileInputStream(mFile);
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
