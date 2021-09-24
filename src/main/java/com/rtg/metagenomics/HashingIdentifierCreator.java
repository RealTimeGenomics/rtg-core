/*
 * Copyright (c) 2017. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */
package com.rtg.metagenomics;

import java.util.Random;

import com.rtg.mode.DnaUtils;

/**
 * Identifier creator using read name and hash of sequence bases.  The hashing is done
 * in such a way the forward and reverse complement versions of the sequence hash to
 * the same value (provided we are told the direction).
 */
public class HashingIdentifierCreator implements IdentifierCreator {

  /** Randomly generated arrays used to compute <code>irvineHash</code> codes */
  private static final long[] HASH_BLOCKS;
  static {
    HASH_BLOCKS = new long[256];
    final Random r = new Random(1L); // use same seed for deterministic behavior
    for (int i = 0; i < 256; ++i) {
      HASH_BLOCKS[i] = r.nextLong();
    }
  }

  static long irvineHash(final byte[] buf) {
    long hash = 0;
    for (int k = 0; k < buf.length; ++k) {
      hash = Long.rotateLeft(hash, 1) ^ HASH_BLOCKS[(buf[k] + k) & 0xFF];
    }
    return hash;
  }

  static long irvineHashRC(final byte[] buf) {
    long hash = 0;
    for (int k = 0, j = buf.length - 1; j >= 0; ++k, --j) {
      hash = Long.rotateLeft(hash, 1) ^ HASH_BLOCKS[(DnaUtils.complement(buf[j]) + k) & 0xFF];
    }
    return hash;
  }

  @Override
  public long getIdentifier(final String readName, final byte[] readBases, final boolean rc) {
    return readName.hashCode() ^ (rc ? irvineHashRC(readBases) : irvineHash(readBases));
  }
}
