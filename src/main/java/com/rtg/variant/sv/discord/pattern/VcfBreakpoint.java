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

package com.rtg.variant.sv.discord.pattern;

import com.rtg.util.CompareHelper;
import com.rtg.util.Utils;
import com.rtg.vcf.BreakpointAlt;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.VcfUtils;

/**
 * After reading breakpoints in from a VCF file abstract out the information we actually care about
 */
public class VcfBreakpoint implements Comparable<VcfBreakpoint> {

  private final String mLocalChr;
  private final int mLocalPos;

  private final BreakpointAlt mBreakpointAlt;

  private final int mDepth;

  /**
   * Create a breakpoint from a VCF record, assuming that the first ALT contains the breakend.
   * @param record the VCF record that describes a breakpoint
   */
  public VcfBreakpoint(VcfRecord record) {
    this(record, new BreakpointAlt(record.getAltCalls().get(0)));
  }

  private VcfBreakpoint(VcfRecord record, BreakpointAlt alt) {
    mLocalChr = record.getSequenceName();
    mLocalPos = record.getStart();
    mBreakpointAlt = alt;
    mDepth = countDepth(record);
    assert record.getRefCall().length() == 1 : "We only handle breakends with single-base reference span";
    assert record.getRefCall().length() == alt.getRefSubstitution().length() : "We only handle breakends with same substitution length as reference span";
  }

  /**
   * Create a breakpoint without a VCF record
   * @param localChr which reference sequence is the local end of the breakpoint in
   * @param localPos position within the reference sequence of the local end
   * @param remoteChr which reference sequence is the remote end of the breakpoint in
   * @param remotePos position within the reference sequence of the remote end
   * @param localUp which direction is the local breakpoint
   * @param remoteUp which direction is the remote breakpoint
   * @param depth number of reads spanning the breakpoint
   */
  public VcfBreakpoint(String localChr, int localPos, String remoteChr, int remotePos, boolean localUp, boolean remoteUp, int depth) {
    mLocalChr = localChr;
    mLocalPos = localPos;
    mBreakpointAlt = new BreakpointAlt("N", localUp, remoteChr, remotePos, remoteUp);
    mDepth = depth;
  }

  private static int countDepth(VcfRecord record) {
    final String depth = record.getInfo(VcfUtils.INFO_COMBINED_DEPTH);
    return depth == null ? 0 : Integer.parseInt(depth);
  }

  public String getLocalChr() {
    return mLocalChr;
  }
  public int getLocalPos() {
    return mLocalPos;
  }
  public boolean isLocalUp() {
    return mBreakpointAlt.isLocalUp();
  }
  public boolean isRemoteUp() {
    return mBreakpointAlt.isRemoteUp();
  }
  public int getRemotePos() {
    return mBreakpointAlt.getRemotePos();
  }
  public String getRemoteChr() {
    return mBreakpointAlt.getRemoteChr();
  }

  /**
   * @return the number of reads contributing to this breakpoint
   */
  public int getDepth() {
    return mDepth;
  }

  @Override
  public int compareTo(VcfBreakpoint o) {
    return new CompareHelper()
        .compare(mLocalChr, o.mLocalChr)
        .compare(mLocalPos, o.mLocalPos)
        .compare(getRemoteChr(), o.getRemoteChr())
        .compare(getRemotePos(), o.getRemotePos())
        .compare(isLocalUp(), o.isLocalUp())
        .compare(isRemoteUp(), o.isRemoteUp())
        .result();
  }

  @Override
  public boolean equals(Object obj) {
    return obj != null && obj.getClass().equals(this.getClass()) && this.compareTo((VcfBreakpoint) obj) == 0;
  }

  @Override
  public int hashCode() {
    return Utils.hash(new Object[]{mLocalChr, mLocalPos, getRemoteChr(), getRemotePos(), isLocalUp(), isRemoteUp()});
  }

  @Override
  public String toString() {
    return "VcfBreakpoint: " + mLocalChr + " " + mLocalPos + " " + mBreakpointAlt;
  }

}
