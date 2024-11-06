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
package com.rtg.index.hash;

import com.rtg.util.integrity.Exam;

/**
 */
public final class PrimeUtils {

  private PrimeUtils() { } //to prevent instantiation

  /**
   * Get a prime number suitable for use as a hash multiplier for winds of length bits.
   * Each is &lt; 2^b a prime and ~ 2/3 * 2^b (it is not clear if this last is necessary or even helps).
   * @param bits number of bits in window.
   * @return suitable multiplier for hash function.
   */
  public static long prime(final int bits) {
    if (bits < 1 || bits > 64) {
      throw new IllegalArgumentException("Number of bits should be between 1 and 64 inclusive. bits=" + bits);
    }
    return PRIME[bits];
  }

  /**
   * Primes for each bit length.
   * Each is &lt; 2^b a prime and ~ 2/3 * 2^b (it is not clear if this last is necessary or even helps)
   */
  private static final long[] PRIME = new long[65];
  /** Inverses of the primes used for multiplying. */
  private static final long[] INVERSE = new long[65];
  static {
    PRIME[1] = 1L;
    PRIME[2] = 3L;
    PRIME[3] = 3L;
    PRIME[4] = 11L;
    PRIME[5] = 23L;
    PRIME[6] = 43L;
    PRIME[7] = 89L;
    PRIME[8] = 173L;
    PRIME[9] = 347L;
    PRIME[10] = 683L;
    PRIME[11] = 1367L;
    PRIME[12] = 2731L;
    PRIME[13] = 5471L;
    PRIME[14] = 10937L;
    PRIME[15] = 21851L;
    PRIME[16] = 43691L;
    PRIME[17] = 87383L;
    PRIME[18] = 174763L;
    PRIME[19] = 349529L;
    PRIME[20] = 699053L;
    PRIME[21] = 1398107L;
    PRIME[22] = 2796203L;
    PRIME[23] = 5592407L;
    PRIME[24] = 11184829L;
    PRIME[25] = 22369661L;
    PRIME[26] = 44739259L;
    PRIME[27] = 89478503L;
    PRIME[28] = 178956983L;
    PRIME[29] = 357913951L;
    PRIME[30] = 715827883L;
    PRIME[31] = 1431655777L;
    PRIME[32] = 2863311551L;
    PRIME[33] = 5726623081L;
    PRIME[34] = 11453246153L;
    PRIME[35] = 22906492261L;
    PRIME[36] = 45812984533L;
    PRIME[37] = 91625969003L;
    PRIME[38] = 183251937977L;
    PRIME[39] = 366503875943L;
    PRIME[40] = 733007751857L;
    PRIME[41] = 1466015503703L;
    PRIME[42] = 2932031007403L;
    PRIME[43] = 5864062014871L;
    PRIME[44] = 11728124029663L;
    PRIME[45] = 23456248059227L;
    PRIME[46] = 46912496118533L;
    PRIME[47] = 93824992237007L;
    PRIME[48] = 187649984473789L;
    PRIME[49] = 375299968947551L;
    PRIME[50] = 750599937895091L;
    PRIME[51] = 1501199875790171L;
    PRIME[52] = 3002399751580363L;
    PRIME[53] = 6004799503160669L;
    PRIME[54] = 12009599006321351L;
    PRIME[55] = 24019198012642661L;
    PRIME[56] = 48038396025285301L;
    PRIME[57] = 96076792050570629L;
    PRIME[58] = 192153584101141199L;
    PRIME[59] = 384307168202282401L;
    PRIME[60] = 768614336404564651L;
    PRIME[61] = 1537228672809129329L;
    PRIME[62] = 3074457345618258637L;
    PRIME[63] = 6148914691236517223L;
    PRIME[64] = 6148914691236517223L;
    INVERSE[1] = 1L;
    INVERSE[2] = 3L;
    INVERSE[3] = 3L;
    INVERSE[4] = 3L;
    INVERSE[5] = 7L;
    INVERSE[6] = 3L;
    INVERSE[7] = 105L;
    INVERSE[8] = 37L;
    INVERSE[9] = 211L;
    INVERSE[10] = 3L;
    INVERSE[11] = 1639L;
    INVERSE[12] = 3L;
    INVERSE[13] = 1695L;
    INVERSE[14] = 15241L;
    INVERSE[15] = 21203L;
    INVERSE[16] = 3L;
    INVERSE[17] = 26215L;
    INVERSE[18] = 3L;
    INVERSE[19] = 238313L;
    INVERSE[20] = 149797L;
    INVERSE[21] = 1233619L;
    INVERSE[22] = 3L;
    INVERSE[23] = 6710887L;
    INVERSE[24] = 5185685L;
    INVERSE[25] = 16354261L;
    INVERSE[26] = 39717491L;
    INVERSE[27] = 70907479L;
    INVERSE[28] = 65295111L;
    INVERSE[29] = 240666271L;
    INVERSE[30] = 3L;
    INVERSE[31] = 552210081L;
    INVERSE[32] = 1126548799L;
    INVERSE[33] = 3785394905L;
    INVERSE[34] = 14914391929L;
    INVERSE[35] = 9503757421L;
    INVERSE[36] = 33548091005L;
    INVERSE[37] = 65547808579L;
    INVERSE[38] = 242915359625L;
    INVERSE[39] = 165964019287L;
    INVERSE[40] = 57869033041L;
    INVERSE[41] = 439804651111L;
    INVERSE[42] = 3L;
    INVERSE[43] = 3125515287079L;
    INVERSE[44] = 1904886386975L;
    INVERSE[45] = 20696689464019L;
    INVERSE[46] = 46479724014029L;
    INVERSE[47] = 20821436633391L;
    INVERSE[48] = 163767259177109L;
    INVERSE[49] = 427065481905823L;
    INVERSE[50] = 135107988821115L;
    INVERSE[51] = 794752875418323L;
    INVERSE[52] = 835719518481123L;
    INVERSE[53] = 4307790947919605L;
    INVERSE[54] = 15471189308143351L;
    INVERSE[55] = 7665701493396589L;
    INVERSE[56] = 32542139242935197L;
    INVERSE[57] = 119928023643544397L;
    INVERSE[58] = 18510207642770479L;
    INVERSE[59] = 86342139111525985L;
    INVERSE[60] = 3L;
    INVERSE[61] = 750093509021322129L;
    INVERSE[62] = 1164114917272932869L;
    INVERSE[63] = 7309087274488690263L;
    INVERSE[64] = 7309087274488690263L;
  }

  /**
   * Get a prime number suitable for use as a hash multiplier for winds of length bits.
   * Each is &lt; 2^b a prime and ~ 2/3 * 2^b (it is not clear if this last is necessary or even helps).
   * @param bits number of bits in window.
   * @return suitable multiplier for hash function.
   */
  public static long primeInverse(final int bits) {
    if (bits < 1 || bits > 64) {
      throw new IllegalArgumentException("Number of bits should be between 1 and 64 inclusive. bits=" + bits);
    }
    return INVERSE[bits];
  }

  static void integrity() {
    Exam.assertTrue(PRIME.length == 65);
    Exam.assertTrue(INVERSE.length == 65);
  }
}

