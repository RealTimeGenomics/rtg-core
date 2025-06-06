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

package com.rtg.variant.sv.discord;

import java.io.Serializable;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.List;

import com.rtg.util.StringUtils;
import com.rtg.variant.sv.bndeval.AbstractBreakpointGeometry;

import htsjdk.samtools.SAMRecord;

/**
 */
class DiscordantReadSet {

  private final String mSequenceName;
  private final int mMaxVariation;
  private final List<BreakpointConstraint> mConstraints = new LinkedList<>();
  private final LinkedList<SAMRecord> mRecords = new LinkedList<>();
  private BreakpointConstraint mIntersection;
  private BreakpointConstraint mUnion;

  /**
   * @param referenceName name of sequence.
   * @param maxVariation maximum distance that can occur between current position and an overlapping <code>DiscordantReadSet</code>
   * @param record first constraint to be added to the set (necessary to avoid all sorts of special cases with empty set).
   */
  DiscordantReadSet(String referenceName, int maxVariation, BreakpointConstraint record) {
    mSequenceName = referenceName;
    mMaxVariation = maxVariation;
    mIntersection = record;
    mUnion = record;
    mConstraints.add(record);
  }

  /**
   * @return the earliest position along this sequence covered by this read set.
   */
  int unionPosition() {
    return mUnion.position(mSequenceName);
  }

  /**
   * @return start position for a constraint beyond which nothing can overlap this read set.
   */
  int flushPosition() {
    return mUnion.getXLo() + mMaxVariation;
  }


  /**
   * Check if the geometry overlaps with this read set.
   * @param bg geometry to be checked.
   * @return true iff the geometry overlaps this read set.
   */
  boolean belongs(AbstractBreakpointGeometry bg) {
    if (!mUnion.overlap(bg)) {
      return false;
    }
    if (bg.getXName().equals(bg.getYName())) {
      final boolean me = mUnion.getXLo() <= mUnion.getYLo();
      final boolean you = bg.getXLo() <= bg.getYLo();
      if (me != you) {
        return false;
      }
    }
    if (mIntersection != null && mIntersection.overlap(bg)) {
      return true;
    }
    for (final AbstractBreakpointGeometry bg1 : mConstraints) {
      if (bg.overlap(bg1)) {
        return true;
      }
    }
    return false;
  }

  /**
   * Add the geometry to this read set.
   * @param bg geometry to be added.
   */
  void add(BreakpointConstraint bg) {
    if (mIntersection != null) {
      mIntersection = mIntersection.intersect(bg);
    }
    //final BreakpointConstraint union = mUnion;
    mUnion = mUnion.union(bg);
    //assert union != null : bg + ":" + union + ":" + mUnion;
    mConstraints.add(bg);
  }

  /**
   * Add all members of drs to this read set.
   * @param drs set of geometries to be added.
   */
  void addAll(DiscordantReadSet drs) {
    if (mIntersection != null) {
      if (drs.mIntersection == null) {
        mIntersection = null;
      } else {
        mIntersection = mIntersection.intersect(drs.mIntersection);
      }
    }
    mUnion = mUnion.union(drs.mUnion);
    assert mUnion != null;
    mConstraints.addAll(drs.mConstraints);
    if (mRecords != null) {
      mRecords.addAll(drs.mRecords);
    }
  }

  @Override
  public String toString() {
    final StringBuilder sb = new StringBuilder();
    sb.append("DiscordantReadSet:").append(StringUtils.LS);
    sb.append("union=").append(mUnion).append(StringUtils.LS);
    sb.append("intersection=").append(mIntersection).append(StringUtils.LS);
    sb.append("constraints count=").append(mConstraints.size());
    sb.append(StringUtils.LS);
    for (final AbstractBreakpointGeometry bg : mConstraints) {
      sb.append("    ").append(bg).append(StringUtils.LS);
    }
    return sb.toString();
  }

  /**
   * @return number of members of this read set.
   */
  int getCounts() {
    return mConstraints.size();
  }

  /**
   * @return a minimal geometry that encloses all the members of this read set.
   */
  BreakpointConstraint getUnion() {
    return mUnion;
  }

  /**
   * @return an intersection of all the members of this read set (null if the intersection is empty).
   */
  BreakpointConstraint getIntersection() {
    return mIntersection;
  }

  /**
   * @return the name of the current sequence (the one used to construct this read set).
   */
  String getSequenceName() {
    return mSequenceName;
  }

  static class FlushPositionComparator implements Comparator<DiscordantReadSet>, Serializable {
    @Override
    public int compare(DiscordantReadSet arg0, DiscordantReadSet arg1) {
      final int p0 = arg0.unionPosition();
      final int p1 = arg1.unionPosition();
      // Important - we only ever want to return 0 for the same exact object, as the results are being stored in a sorted set
      return p0 == p1 ? Integer.compare(System.identityHashCode(arg0), System.identityHashCode(arg1)) : Integer.compare(p0, p1);
    }
  }

  void add(SAMRecord record) {
    if (mRecords != null) {
      mRecords.add(record);
    }
  }

  List<SAMRecord> getRecords() {
    return mRecords;
  }

}
