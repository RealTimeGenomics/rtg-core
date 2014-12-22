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

package com.rtg.variant.sv.discord;

import java.io.Serializable;
import java.util.Comparator;

import com.rtg.util.CompareHelper;

/**
*         Date: 1/03/12
*         Time: 4:02 PM
*/
class DebugComparator implements Comparator<DiscordantReadSet>, Serializable {
  @Override
  public int compare(DiscordantReadSet o1, DiscordantReadSet o2) {
    return new CompareHelper()
        .compare(o1.getUnion().getX(), o2.getUnion().getX())
        .compare(o1.getUnion().getZ(), o2.getUnion().getZ())
        .compare(o1.getUnion().getYName(), o2.getUnion().getYName())
        .compare(o1.getUnion().getY(), o2.getUnion().getY())
        .compare(o1.getUnion().getW(), o2.getUnion().getW())
        .compare(o1.getCounts(), o2.getCounts())
        .result();
  }
}
