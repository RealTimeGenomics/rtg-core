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
package com.rtg.index.params;

import com.rtg.util.StringUtils;


/**
 * Utilities of general use to the <code>Params</code> classes.
 */
public final class ParamsUtils {

  private ParamsUtils() { } //prevent instantiation

  /**
   * Create a single line with a memory display (used by all the parameter classes).
   * @param id identifier for the memory being displayed (must not contain and white space).
   * @param bytes number of bytes consumed by this item of memory.
   * @param params optional parameters for specifying sizes of things.
   * @return a single line with a memory display.
   */
  public static String memToString(final String id, final long bytes, final long...params) {
    if (id.matches(".*\\s.*")) {
      throw new RuntimeException();
    }
    final StringBuilder sb = new StringBuilder();
    sb.append("\tMemory\t").append(id).append("\t").append(StringUtils.commas(bytes));
    for (long p : params) {
      sb.append("\t").append(StringUtils.commas(p));
    }
    sb.append(StringUtils.LS);
    return sb.toString();
  }

}

