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

import java.io.IOException;

/**
 * Just like <code>java.lang.Runnable</code> except that we are running
 * I/O heavy code, so expect that <code>java.io.IOException</code> may
 * be thrown.
 *
 */
public interface IORunnable {

  /**
   * Like <code>java.io.Runnable.run()</code>, but it is allowed to
   * throw IOException.
   *
   * @exception IOException if an error occurs.
   */
  void run() throws IOException;
}
