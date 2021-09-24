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

package com.rtg.assembler.graph;

/**
 * Represent a contig as an array of bytes
 */
public class ContigByte implements Contig {
  byte[] mBytes;

  /**
   * Construct a contig containing the given bytes
   * @param bytes the bytes representing contig bases
   */
  public ContigByte(byte[] bytes) {
    mBytes = bytes;
  }

  @Override
  public int length() {
    return mBytes.length;
  }

  @Override
  public byte nt(int index) {
    return mBytes[index];
  }
}
