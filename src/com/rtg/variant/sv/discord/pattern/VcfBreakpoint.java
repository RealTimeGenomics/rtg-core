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

package com.rtg.variant.sv.discord.pattern;

import java.util.Collection;

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
    int depth = 0;
    final Collection<String> depths = record.getInfo().get(VcfUtils.INFO_COMBINED_DEPTH);
    if (depths != null) {
      for (String str : depths) {
        depth += Integer.parseInt(str);
      }
    }
    return depth;
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
