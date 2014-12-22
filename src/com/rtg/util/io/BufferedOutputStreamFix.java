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
package com.rtg.util.io;

import java.io.BufferedOutputStream;
import java.io.IOException;
import java.io.OutputStream;

/**
 * {@link BufferedOutputStream} extends {@link java.io.FilterOutputStream} which has a serious bug in it's close
 * method (silently ignores IOException). This class provides a working version.
 */
public class BufferedOutputStreamFix extends BufferedOutputStream {

  private boolean mClosed = false;

  /**
   * @param out {@link BufferedOutputStream#BufferedOutputStream(java.io.OutputStream)}
   */
  public BufferedOutputStreamFix(OutputStream out) {
    super(out);
  }

  /**
   * @param out {@link BufferedOutputStream#BufferedOutputStream(java.io.OutputStream)}
   * @param size {@link BufferedOutputStream#BufferedOutputStream(java.io.OutputStream)}
   */
  public BufferedOutputStreamFix(OutputStream out, int size) {
    super(out, size);
  }

  @Override
  public void close() throws IOException {
    if (!mClosed) {
      mClosed = true;
      try {
        flush();
      } finally {
        out.close();
      }
    }
  }
}
