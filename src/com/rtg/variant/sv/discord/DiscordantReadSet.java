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
import java.util.LinkedList;
import java.util.List;

import com.rtg.util.StringUtils;

import htsjdk.samtools.SAMRecord;

/**
 */
class DiscordantReadSet {


  final String mSequenceName;
  private final int mMaxVariation;
  private BreakpointConstraint mIntersection;
  private BreakpointConstraint mUnion;
  final List<BreakpointConstraint> mConstraints = new LinkedList<>();

  final LinkedList<SAMRecord> mRecords = !DiscordantTool.DUMP_RECORDS.equals(DiscordantTool.NONE) ? new LinkedList<SAMRecord>() : null;

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
    return mUnion.getX() + mMaxVariation;
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
      final boolean me = mUnion.getX() <= mUnion.getY();
      final boolean you = bg.getX() <= bg.getY();
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
    //TODO make more efficient
    //    for (final BreakpointConstraint bg : drs.mConstraints) {
    //      add(bg);
    //    }
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
  }

  @Override
  public String toString() {
    final StringBuilder sb = new StringBuilder();
    sb.append("DiscordantReadSet:").append(StringUtils.LS);
    sb.append("union=").append(mUnion).append(StringUtils.LS);
    sb.append("intersection=").append(mIntersection == null ? null : mIntersection).append(StringUtils.LS);
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
    if (record != null && mRecords != null) {
      mRecords.add(record);
    }
  }

  List<SAMRecord> getRecords() {
    return mRecords;
  }

}
