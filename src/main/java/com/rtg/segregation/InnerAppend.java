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

package com.rtg.segregation;

/**
 * Helper for printing set style output.
 */
public class InnerAppend {

  /**
   * @return suitable for outputting a set.
   */
  public static InnerAppend innerSet() {
    return new InnerAppend("{", ", ", "}");
  }

  private final String mPrefix;
  private final String mSeparator;
  private final String mSuffix;

  private boolean mFirstDone = false;

  /**
   * @param prefix done once at start.
   * @param separator before each appended item except the first.
   * @param suffix done once at end.
   */
  public InnerAppend(String prefix, String separator, String suffix) {
    mPrefix = prefix;
    mSeparator = separator;
    mSuffix = suffix;
  }

  /**
   * @return the prefix.
   */
  public String start() {
    assert !mFirstDone;
    return mPrefix;
  }

  /**
   * @param str to be returned with any necessary prepended constant.
   * @return str with prepended constant as necessary.
   */
  public String inner(final String str) {
    if (mFirstDone) {
      return mSeparator + str;
    }
    mFirstDone = true;
    return str;
  }

  /**
   * @return the suffix.
   */
  public String end() {
    return mSuffix;
  }
}
