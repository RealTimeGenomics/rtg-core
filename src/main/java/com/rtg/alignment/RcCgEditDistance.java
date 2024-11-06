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
