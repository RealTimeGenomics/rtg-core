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

import com.reeltwo.jumble.annotations.TestClass;

/**
 * <code>AnonymousFlag</code> is a flag with no name.
 */
@TestClass(value = {"com.rtg.util.cli.CFlagsTest"})
public class AnonymousFlag extends Flag {

  private static int sAnonCounter = 0;
  private static synchronized int nextAnonCounter() {
    return ++sAnonCounter;
  }

  /** This specifies the ordering. */
  private final int mFlagRank;

  /**
   * Constructor for an anonymous <code>Flag</code>. These flags aren't
   * referred to by name on the command line -- their values are assigned
   * based on their position in the command line.
   *
   * @param flagDescription a name used when printing help messages.
   * @param paramType a <code>Class</code> denoting the type of values to be
   * accepted.
   * @param paramDescription a description of the meaning of the flag.
   */
  public AnonymousFlag(final String flagDescription, final Class<?> paramType,
      final String paramDescription) {
    super(null, null, flagDescription, 1, 1, paramType, paramDescription, null, "");
    mFlagRank = nextAnonCounter();
  }

  @Override
  String getFlagUsage() {
    final StringBuilder sb = new StringBuilder();
    sb.append(getParameterDescription());
    if (getMaxCount() > 1) {
      sb.append('+');
    }
    return sb.toString();
  }

  @Override
  String getCompactFlagUsage() {
    return getFlagUsage();
  }

  @Override
  public boolean equals(final Object other) {
    return super.equals(other);
  }

  @Override
  public int hashCode() {
    return super.hashCode();
  }

  @Override
  public int compareTo(final Flag other) {
    if (other instanceof AnonymousFlag) {
      return mFlagRank - ((AnonymousFlag) other).mFlagRank;
    }
    return 1;
  }
}
