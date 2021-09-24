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

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.InputStream;
import java.io.OutputStream;

import com.rtg.util.integrity.IntegralAbstract;

/**
 */
public class StoreFileString extends IntegralAbstract implements StoreFile {

  private final ByteArrayOutputStream mBos = new ByteArrayOutputStream();

  private final String mName;

  /**
   * @param name of this store file.
   */
  public StoreFileString(final String name) {
    mName = name;
  }

  @Override
  public OutputStream outputStream() {
    return mBos;
  }

  @Override
  public InputStream inputStream() {
    return new ByteArrayInputStream(mBos.toByteArray());
  }

  @Override
  public String content() {
    return mBos.toString();
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
