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

package com.rtg.variant.bayes.complex;

import java.util.HashMap;

import com.rtg.variant.realign.AllPaths;
import com.rtg.variant.realign.RealignParams;
import com.rtg.variant.realign.ScoreFastUnderflow;
import com.rtg.variant.realign.ScoreFastUnderflowCG;

/**
 * Keep singleton score interfaces per read group and CG flag.
 */
public final class ScoreInterfaceMemo implements ScoreInterfaceMemoInterface {

  final HashMap<RealignParams, AllPaths> mCache = new HashMap<>();
  final HashMap<RealignParams, AllPaths> mCacheCG = new HashMap<>();

  @Override
  public AllPaths getScoreInterface(final RealignParams params) {
    final AllPaths s;
    if (EvidenceComplex.CG_ALLPATHS && params.machineType() != null && params.machineType().isCG()) {
      s = mCacheCG.computeIfAbsent(params, ScoreFastUnderflowCG::new);
    } else {
      s = mCache.computeIfAbsent(params, ScoreFastUnderflow::new);
    }
    return s;
  }
}
