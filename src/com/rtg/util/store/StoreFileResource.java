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

import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.Reader;

import com.rtg.util.Resources;
import com.rtg.util.integrity.IntegralAbstract;
import com.rtg.util.io.FileUtils;

/**
 */
public class StoreFileResource extends IntegralAbstract implements StoreFile {

  private final String mResource;

  private final String mName;

  /**
   * @param resource name of parent resource.
   * @param name of child of the resource.
   */
  StoreFileResource(final String resource, final String name) {
    mResource = resource;
    mName = name;
  }

  @Override
  public OutputStream outputStream() {
    throw new UnsupportedOperationException();
  }

  @Override
  public InputStream inputStream() {
    final String name = mResource + "/" + mName;
    return Resources.getResourceAsStream(name);
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
