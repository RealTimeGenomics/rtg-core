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
package com.rtg.util.cli;

/**
 * An interface for objects that perform custom flag validation across a set
 * of flags.
 */
public interface Validator {

  /**
   * Returns false if the current flag settings are not valid. An error
   * message should be supplied by calling <code>flags.setParseMessage()</code>.
   *
   * @param flags a <code>CFlags</code>.
   * @return true if the current flag settings are valid.
   */
  boolean isValid(CFlags flags);
}
