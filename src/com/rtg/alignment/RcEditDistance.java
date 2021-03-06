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
package com.rtg.alignment;

import com.reeltwo.jumble.annotations.JumbleIgnore;
import com.rtg.mode.DNA;
import com.rtg.mode.DnaUtils;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * An edit distance wrapper that stores an internal reverse complement version
 * of the template, and always calls the wrapped unidirectional edit distance
 * implementation(s).
 */
public class RcEditDistance implements BidirectionalEditDistance {

  private UnidirectionalEditDistance mEd;
  private UnidirectionalEditDistance mEdRC;
  private byte[] mTemplate = null;
  private byte[] mTemplateRC = null;

  /**
   * Create a new RcEditDistance
   * @param ed the unidirectional edit distance implementation to wrap
   */
  public RcEditDistance(final UnidirectionalEditDistance ed) {
    mEd = ed;
    mEdRC = null;
  }

  /**
   * Create a new RcEditDistance, with two wrapped implementations, one for
   * forward and one for reverse. This is useful in the case that the wrapped
   * implementation caches data based on the template.
   * @param ed the unidirectional edit distance implementation to wrap for the forward case
   * @param edRev the unidirectional edit distance implementation to wrap for the reverse case
   */
  public RcEditDistance(UnidirectionalEditDistance ed, UnidirectionalEditDistance edRev) {
    mEd = ed;
    mEdRC = edRev;
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
    if (mEdRC != null && rc) {
      res = mEdRC.calculateEditDistance(read, rlen, mTemplateRC, mTemplateRC.length - rlen - zeroBasedStart, maxScore, maxShift, cgLeft);
    } else {
      res = mEd.calculateEditDistance(read, rlen, rc ? mTemplateRC : template, rc ? mTemplateRC.length - rlen - zeroBasedStart : zeroBasedStart, maxScore, maxShift, cgLeft);
    }
    //need to translate the position back in reverse complementation
    if (rc && res != null && res[ActionsHelper.ALIGNMENT_SCORE_INDEX] != Integer.MAX_VALUE) {
        res[ActionsHelper.TEMPLATE_START_INDEX] = mTemplateRC.length - ActionsHelper.zeroBasedTemplateEndPos(res);
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
    if (mEdRC != null) {
      Diagnostic.developerLog("RcEditDistance: sub ED RC");
      mEdRC.logStats();
    }
    mTemplate = null;
    mTemplateRC = null;
    mEd = null;
    mEdRC = null;
  }
}
