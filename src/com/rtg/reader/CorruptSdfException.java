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

import java.io.File;
import java.io.IOException;

/**
 * Indicates an error while trying to read SDF
 */
public class CorruptSdfException extends IOException {

  /**
   * @see IOException#IOException(java.lang.String)
   * @param message as in IOException
   */
  public CorruptSdfException(String message) {
    super(message);
  }

  /**
   * @see IOException#IOException(java.lang.String)
   * @param sdfdir the File of the corrupt SDF
   */
  public CorruptSdfException(File sdfdir) {
    super("The SDF " + sdfdir + " is invalid");
  }

  /**
   * @see IOException#IOException()
   */
  public CorruptSdfException() {
    this("Invalid SDF");
  }

}
