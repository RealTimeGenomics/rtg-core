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
import com.rtg.vcf.VcfRecord;

/**
 * After reading breakpoints in from a VCF file abstract out the information we actually care about
 *         Date: 16/03/12
 *         Time: 10:16 AM
 */
public class VcfBreakpoint implements Comparable<VcfBreakpoint> {
  private final String mLocalChr;
  private final int mLocalPos;
  private final String mRemoteChr;
  private final int mRemotePos;
  private final boolean mLocalDirection;
  private final boolean mRemoteDirection;
  private final int mDepth;

  public boolean isLocalUp() {
    return mLocalDirection;
  }
  public boolean isRemoteUp() {
    return mRemoteDirection;
  }

  public int getRemotePos() {
    return mRemotePos;
  }

  public String getLocalChr() {
    return mLocalChr;
  }

  public int getLocalPos() {
    return mLocalPos;
  }

  public String getRemoteChr() {
    return mRemoteChr;
  }

  /**
   * Create a breakpoint from a VCF record
   * @param  record the VCF record that describes a breakpoint
   */
  public VcfBreakpoint(VcfRecord record) {
    this(record, new BreakpointAlt(record.getAltCalls().get(0)));
  }
  static int countDepth(VcfRecord record) {
    int depth = 0;
    final Collection<String> depths = record.getInfo().get("DP");
    if (depths != null) {
      for (String str : depths) {
        depth += Integer.parseInt(str);
      }
    }
    return depth;
  }
  VcfBreakpoint(VcfRecord record, BreakpointAlt alt) {
    this(record.getSequenceName(), record.getOneBasedStart(), alt.getRemoteChr(), alt.getRemotePos(), alt.isLocalUp(), alt.isRemoteUp(), countDepth(record));

  }

  /**
   * Create a breakpoint without a VCF record
   * @param localChr which reference sequence is the local end of the breakpoint in
   * @param localPos position within the reference sequence of the local end, 1-based
   * @param remoteChr which reference sequence is the remote end of the breakpoint in
   * @param remotePos position within the reference sequence of the remote end
   * @param localDirection which direction is the local breakpoint
   * @param remoteDirection which direction is the remote breakpoint
   * @param depth number of reads spanning the breakpoint
   */
  public VcfBreakpoint(String localChr, int localPos, String remoteChr, int remotePos, boolean localDirection, boolean remoteDirection, int depth) {
    mRemoteChr = remoteChr;
    mLocalPos = localPos;
    mRemotePos = remotePos;
    mLocalChr = localChr;
    mLocalDirection = localDirection;
    mRemoteDirection = remoteDirection;
    mDepth = depth;
  }

  @Override
  public int compareTo(VcfBreakpoint o) {
    return new CompareHelper()
        .compare(mLocalChr, o.mLocalChr)
        .compare(mLocalPos, o.mLocalPos)
        .compare(mRemoteChr, o.mRemoteChr)
        .compare(mRemotePos, o.mRemotePos)
        .compare(mLocalDirection, o.mLocalDirection)
        .compare(mRemoteDirection, o.mRemoteDirection)
        .result();
  }

  @Override
  public boolean equals(Object obj) {
    return obj != null && obj.getClass().equals(this.getClass()) && this.compareTo((VcfBreakpoint) obj) == 0;
  }

  @Override
  public int hashCode() {
    return Utils.hash(new Object[]{mLocalChr, mLocalPos, mRemoteChr, mRemotePos, mLocalDirection, mRemoteDirection});
  }

  @Override
  public String toString() {
    return "VcfBreakpoint: " + mLocalChr + " " + mLocalPos + " " + mRemoteChr + " " + mRemotePos + " " + mLocalDirection + " " + mRemoteDirection;
  }

  /**
   * @return the number of reads contributing to this breakpoint
   */
  public int getDepth() {
    return mDepth;
  }
}
