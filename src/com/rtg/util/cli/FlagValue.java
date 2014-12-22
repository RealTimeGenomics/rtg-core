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
 * Encapsulates a flag and value pairing. This is used when retrieving the set
 * of flags in the order they were set.
 */
public class FlagValue {
  private Flag mFlag;

  private Object mValue;

  FlagValue(final Flag flag, final Object value) {
    mFlag = flag;
    mValue = value;
  }

  /**
   * Gets the Flag that this value was supplied to.
   *
   * @return the Flag that this value was supplied to.
   */
  public Flag getFlag() {
    return mFlag;
  }

  /**
   * Gets the value supplied to the flag.
   *
   * @return the value supplied to the flag.
   */
  public Object getValue() {
    return mValue;
  }

  /**
   * Gets a human-readable description of the flag value.
   *
   * @return a human-readable description of the flag value.
   */
  public String toString() {
    String name = mFlag.getName();
    if (name == null) {
      name = mFlag.getParameterDescription();
    }
    return name + "=" + mValue;
  }
}
