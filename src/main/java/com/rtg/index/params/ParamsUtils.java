/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

