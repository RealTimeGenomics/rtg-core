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
import com.rtg.variant.realign.HomoPolymerParams;
import com.rtg.variant.realign.RealignParams;
import com.rtg.variant.realign.ScoreFastUnderflowHomopolymer;

/**
 */
public final class ScoreInterfaceMemoIonT implements ScoreInterfaceMemoInterface {

  private final HomoPolymerParams mHomoParamsSimple;

  private final HomoPolymerParams mHomoParamsLog;

  final HashMap<RealignParams, AllPaths> mCache = new HashMap<>();


  /**
   * @param homoParams calibration of homopolymer transitions used in the all paths matrix (using simple arithmetic).
   * @param homoParamsLog  calibration of homopolymer transitions used in the all paths matrix (using log arithmetic).
   */
  ScoreInterfaceMemoIonT(HomoPolymerParams homoParams, HomoPolymerParams homoParamsLog) {
    mHomoParamsSimple = homoParams;
    mHomoParamsLog = homoParamsLog;
  }


  @Override
  public AllPaths getScoreInterface(final RealignParams me) {
    if (me.machineType() != null && me.machineType().isCG()) {
      throw new UnsupportedOperationException();
    }
    AllPaths s;
    s = mCache.get(me);
    if (s == null) {
      s = new ScoreFastUnderflowHomopolymer(me, mHomoParamsSimple, mHomoParamsLog);
      mCache.put(me, s);
    }
    return s;
  }
}
