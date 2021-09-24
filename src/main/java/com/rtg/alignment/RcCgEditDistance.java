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
package com.rtg.alignment;

import com.reeltwo.jumble.annotations.JumbleIgnore;
import com.rtg.mode.DNA;
import com.rtg.mode.DnaUtils;
import com.rtg.reader.CgUtils;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * An edit distance wrapper that stores an internal reverse complement version
 * of the template, and always calls the wrapped unidirectional edit distance
 * implementation(s).
 */
public class RcCgEditDistance implements BidirectionalEditDistance {

  private UnidirectionalEditDistance mEd;
  private byte[] mTemplate = null;
  private byte[] mTemplateRC = null;

  /**
   * Create a new RcEditDistance
   * @param ed the unidirectional edit distance implementation to wrap
   */
  public RcCgEditDistance(final CgGotohEditDistance ed) {
    mEd = ed;
  }

  @Override
  public int[] calculateEditDistance(byte[] read, int rlen, byte[] template, int zeroBasedStart, boolean rc, int maxScore, int maxShift, boolean cgLeft) {
    if (rlen < 0 || zeroBasedStart >= template.length) {
      Diagnostic.developerLog("PROBLEM: RcEditDistance problem: " + rlen + " " + zeroBasedStart + "/" + template.length + " " + maxScore + " " + DnaUtils.bytesToSequenceIncCG(read));
    }
    if (template != mTemplate) {
      mTemplate = template;
      mTemplateRC = template.clone();
      DNA.reverseComplementInPlace(mTemplateRC, 0, mTemplateRC.length);
    }
    final int[] res;
    // We have to correct for the default gap in CG v1 reads.
    final int nastyCGHack = rlen == CgUtils.CG_RAW_READ_LENGTH ? CgGotohEditDistance.CG_INSERT_REGION_DEFAULT_SIZE : 0;
    res = mEd.calculateEditDistance(read, rlen, rc ? mTemplateRC : template, rc ? mTemplateRC.length - rlen - zeroBasedStart - nastyCGHack : zeroBasedStart, maxScore, maxShift, cgLeft);
    //need to translate the position back in reverse complementation
    if (res != null && rc) {
      res[ActionsHelper.TEMPLATE_START_INDEX] = mTemplateRC.length - res[ActionsHelper.TEMPLATE_START_INDEX] - (rlen + ActionsHelper.indelLength(res));
    }
    return res;
  }

  @Override
  @JumbleIgnore
  public void logStats() {
    if (mEd != null) {
      Diagnostic.developerLog("RcEditDistance: sub ED forward");
      mEd.logStats();
    }
    mTemplate = null;
    mTemplateRC = null;
    mEd = null;
  }
}
