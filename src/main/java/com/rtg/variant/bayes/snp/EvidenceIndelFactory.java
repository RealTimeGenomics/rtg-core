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

package com.rtg.variant.bayes.snp;

import com.rtg.launcher.globals.GlobalFlags;
import com.rtg.launcher.globals.CoreGlobalFlags;
import com.rtg.variant.bayes.EvidenceInterface;
import com.rtg.variant.bayes.Model;

/**
 * Factory for indel evidence objects.
 */
public final class EvidenceIndelFactory implements CachedEvidenceFactory {

  //TODO cache these if important
  private static final boolean COMPLEX_REGION_INDEL_EXTENSION = GlobalFlags.isSet(CoreGlobalFlags.COMPLEX_REGION_INDEL_EXTENSION);

  /** Unique instance of factory. */
  public static final CachedEvidenceFactory SINGLETON = new EvidenceIndelFactory();


  // For efficiency precompute the possible evidence objects
  private final EvidenceIndel[] mGoodEvidence;
  private final EvidenceIndel[] mBadEvidence;


  private EvidenceIndelFactory() {
    mGoodEvidence = new EvidenceIndel[EvidenceIndel.DISTINCT_READS];
    mBadEvidence = new EvidenceIndel[EvidenceIndel.DISTINCT_READS];
    for (int i = 0; i < EvidenceIndel.DISTINCT_READS; ++i) {
      mGoodEvidence[i] = new EvidenceIndel(0, i, 0);
      mBadEvidence[i] = new EvidenceIndel(Model.AMBIGUITY_THRESHOLD, i, 0);
    }
  }

  @Override
  public EvidenceInterface evidence(int operationType, int readBasesLeft, int readBasesRight, int mapQ, int phred, int stateIndex, int operationLength, boolean isUnmapped) {
    assert !isUnmapped; // Unmapped evidence should never be coming as indels
    if (mapQ <= Model.AMBIGUITY_PHRED) {
      return COMPLEX_REGION_INDEL_EXTENSION ? new EvidenceIndel(Model.AMBIGUITY_THRESHOLD, operationType, operationLength) : mBadEvidence[operationType];
    } else {
      return COMPLEX_REGION_INDEL_EXTENSION ? new EvidenceIndel(0, operationType, operationLength) : mGoodEvidence[operationType];
    }
  }

  // This one does not care about any of these states.
  @Override
  public int getStateIndex(boolean isForward, boolean isReadPaired, boolean isFirst, boolean isMated) {
    return 0;
  }
}
