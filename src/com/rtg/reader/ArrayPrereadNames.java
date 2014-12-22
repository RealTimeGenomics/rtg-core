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

import java.io.IOException;
import java.io.OutputStream;

/**
 * Utility class for writing names
 */
class ArrayPrereadNames extends PrereadNames {

  private final String[] mNames;

  protected ArrayPrereadNames(final String[] names) {
    mNames = names;
  }

  @Override
  public long length() {
    return mNames.length;
  }

  @Override
  public String name(final long id) {
    return mNames[(int) id];
  }

  @Override
  public void writeName(final Appendable a, final long id) throws IOException {
    a.append(name(id));
  }

  @Override
  public void writeName(final OutputStream os, final long id) throws IOException {
    os.write(name(id).getBytes());
  }
}
