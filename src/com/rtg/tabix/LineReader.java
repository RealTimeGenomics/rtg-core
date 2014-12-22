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
package com.rtg.tabix;

import java.io.Closeable;
import java.io.IOException;

/**
 * Interface for something that can read a line of input.
 *
 */
public interface LineReader extends Closeable {

  /**
   * Reads a line of input
   * @return the read line
   * @throws IOException if a problem occurs during reading
   */
  String readLine() throws IOException;
}
