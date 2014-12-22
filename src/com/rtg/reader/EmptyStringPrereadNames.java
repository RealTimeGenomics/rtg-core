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

import java.io.OutputStream;

/**
 * Names implementation where all values are the empty string
 */
class EmptyStringPrereadNames implements PrereadNamesInterface {
  private final long mLength;

  public EmptyStringPrereadNames(long length) {
    this.mLength = length;
  }

  @Override
  public long length() {
    return mLength;
  }

  @Override
  public String name(long id) {
    return "";
  }

  @Override
  public long calcChecksum() {
    final PrereadHashFunction namef = new PrereadHashFunction();
    for (int k = 0; k < length(); k++) {
      namef.irvineHash("");
      namef.irvineHash(0);
    }
    return namef.getHash();
  }

  @Override
  public long bytes() {
    return 0;
  }

  @Override
  public void writeName(Appendable a, long id) {
  }

  @Override
  public void writeName(OutputStream os, long id) {
  }

}
