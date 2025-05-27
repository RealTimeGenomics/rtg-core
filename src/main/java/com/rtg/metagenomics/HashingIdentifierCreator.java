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
