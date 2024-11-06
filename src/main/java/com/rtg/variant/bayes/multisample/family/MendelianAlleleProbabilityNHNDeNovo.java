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

package com.rtg.variant.bayes.multisample.family;

import com.rtg.variant.bayes.Code;

/**
 * Cases for the Y chromosome in mammals.
 * Mother ignored, father haploid, daughter ignored.
 */
public final class MendelianAlleleProbabilityNHNDeNovo extends MendelianAlleleProbabilityDeNovo {

  /** Father really is the father. */
  public static final MendelianAlleleProbability SINGLETON_HN = new MendelianAlleleProbabilityNHNDeNovo();

  /** Father really is the mother. */
  public static final MendelianAlleleProbability SINGLETON_NH = new MendelianAlleleProbabilityDeNovo() {
    @Override
    public double probabilityLn(Code code, int father, int mother, int child) {
      return SINGLETON_HN.probabilityLn(code, mother, father, child);
    }
  };

  @Override
  public double probabilityLn(Code code, int father, int mother, int child) {
    assert code.homozygous(father);
//    assert mother == -1; //Annonyingly not the case when using FB
//    assert child == -1;
    return Double.NEGATIVE_INFINITY;
  }

}
