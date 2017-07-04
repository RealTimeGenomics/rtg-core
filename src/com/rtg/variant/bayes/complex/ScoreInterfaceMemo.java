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
    if (params.machineType() != null && params.machineType().isCG() && EvidenceComplex.CG_ALLPATHS) {
      s = mCacheCG.computeIfAbsent(params, ScoreFastUnderflowCG::new);
    } else {
      s = mCache.computeIfAbsent(params, ScoreFastUnderflow::new);
    }
    return s;
  }
}
