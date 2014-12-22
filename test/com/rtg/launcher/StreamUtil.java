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
package com.rtg.launcher;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;

/**
 * Handy way of capturing and testing output into a <code>PrintStream</code>
 */
public class StreamUtil extends PrintStream {

  /**
   * Creates the print stream.
   * @return the print stream
   */
  public static PrintStream make() {
    return new StreamUtil(new ByteArrayOutputStream());
  }

  private final ByteArrayOutputStream mOut;

  StreamUtil(final ByteArrayOutputStream baout) {
    super(baout);
    mOut = baout;
  }

  @Override
  public String toString() {
    close();
    try {
      mOut.close();
    } catch (final IOException e) {
      throw new RuntimeException(e);
    }
    return new String(mOut.toByteArray());
  }

}
