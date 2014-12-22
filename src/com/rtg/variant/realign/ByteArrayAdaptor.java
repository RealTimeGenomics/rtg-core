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

package com.rtg.variant.realign;

/**
 */
public abstract class ByteArrayAdaptor {

  /**
   * @param index 0 based.
   * @return the byte at the specified index.
   */
  public abstract byte get(int index);

  /**
   * @return the length of the adapted array.
   */
  public abstract int length();

  @Override
  public String toString() {
    final StringBuilder sb = new StringBuilder();
    sb.append("[");
    for (int i = 0; i < length(); i++) {
      if (i > 0) {
        sb.append(", ");
      }
      sb.append("").append(get(i));
    }
    sb.append("]");
    return sb.toString();
  }


}
