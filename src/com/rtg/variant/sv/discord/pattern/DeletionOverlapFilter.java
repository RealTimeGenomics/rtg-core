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

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import com.rtg.util.AutoAddMap;
import com.rtg.util.CompareHelper;
import com.rtg.variant.sv.discord.DiscordBedRecord;

/**
 * Takes a set of deletion records and returns a non-overlapping subset.
 */
public class DeletionOverlapFilter implements Iterable<DiscordBedRecord> {

  private static class BedRecordAutoAdd extends AutoAddMap<String, ArrayList<DiscordBedRecord>> {
    @Override
    public ArrayList<DiscordBedRecord> make() {
      return new ArrayList<>();
    }
  }

  BedRecordAutoAdd mChrMap = new BedRecordAutoAdd();

  /**
   * Add a deletion
   * @param bed the deletion bed record
   */
  public void add(DiscordBedRecord bed) {
    final ArrayList<DiscordBedRecord> chr = mChrMap.getOrAdd(bed.getSequenceName());
    chr.add(bed);
  }

  @Override
  public Iterator<DiscordBedRecord> iterator() {
    return nonOverlapping().iterator();
  }

  /**
   * Records A &amp; B overlap if A.start &lt; B.end &amp;&amp; A.end &gt; B.start
   * We find a set of records with valid start and union them with valid ends.
   * @return the set of non overlapping deletions
   */
  public List<DiscordBedRecord> nonOverlapping() {
    final List<DiscordBedRecord> result = new ArrayList<>();
    for (final ArrayList<DiscordBedRecord> chrRecords : mChrMap.values()) {
      final TreeSet<DiscordBedRecord> size = new TreeSet<>(new SizeComparator());
      final TreeSet<DiscordBedRecord> start = new TreeSet<>(new StartComparator());
      final TreeSet<DiscordBedRecord> end = new TreeSet<>(new EndComparator());
      size.addAll(chrRecords);
      start.addAll(chrRecords);
      end.addAll(chrRecords);
      while (size.size() > 0) {
        final DiscordBedRecord current = size.pollFirst();

        final Set<DiscordBedRecord> tmp = new HashSet<>();
        tmp.addAll(startSet(start, current));
        tmp.retainAll(endSet(end, current));
        size.removeAll(tmp);
        start.removeAll(tmp);
        end.removeAll(tmp);
        result.add(current);
      }
    }
    return result;
  }

  /**
   * Find a set of records that start before the end of <code>end</code>
   * @param startRecords a set sorted by Bed record start position.
   * @param end find records that start before this ones end.
   * @return all records in <code>startRecords</code> that start before <code>end.getEnd()</code>
   */
  final SortedSet<DiscordBedRecord> startSet(TreeSet<DiscordBedRecord> startRecords, DiscordBedRecord end) {
    final DiscordBedRecord bed = new DiscordBedRecord(end.getSequenceName(), end.getEnd() - 1, Integer.MAX_VALUE);
    return startRecords.headSet(bed);
  }

  /**
   * Find a set of records that end before the start of <code>start</code>
   * @param endRecords a set sorted by Bed record end position.
   * @param start find records that end before this ones start.
   * @return all records in <code>endRecords</code> that start before <code>end.getStart()</code>
   */
  final SortedSet<DiscordBedRecord> endSet(TreeSet<DiscordBedRecord> endRecords, DiscordBedRecord start) {
    final DiscordBedRecord bed = new DiscordBedRecord(start.getSequenceName(), Integer.MIN_VALUE, start.getStart() + 1);
    return endRecords.tailSet(bed);
  }
  class BaseComparator implements Comparator<DiscordBedRecord> {
    @Override
    public int compare(DiscordBedRecord o1, DiscordBedRecord o2) {
      return new CompareHelper()
          .compare(o1.getStart(), o2.getStart())
          .compare(o1.getEnd(), o2.getEnd())
          .compareArray(o1.getAnnotations(), o2.getAnnotations())
          .result();
    }
  }

  class SizeComparator extends BaseComparator {
    @Override
    public int compare(DiscordBedRecord o1, DiscordBedRecord o2) {
      final int res = Integer.compare(Math.abs(o1.getEnd() - o1.getStart()), Math.abs(o2.getEnd() - o2.getStart()));
      if (res != 0) {
        return res;
      }
      return super.compare(o1, o2);
    }
  }
  class StartComparator extends BaseComparator {
    @Override
    public int compare(DiscordBedRecord o1, DiscordBedRecord o2) {
      final int start1 = Math.min(o1.getStart(), o1.getEnd());
      final int start2 = Math.min(o2.getStart(), o2.getEnd());
      final int res = Integer.compare(start1, start2);
      if (res != 0) {
        return res;
      }
      return super.compare(o1, o2);
    }
  }
  class EndComparator extends BaseComparator {
    @Override
    public int compare(DiscordBedRecord o1, DiscordBedRecord o2) {
      final int end1 = Math.max(o1.getStart(), o1.getEnd());
      final int end2 = Math.max(o2.getStart(), o2.getEnd());
      final int res = Integer.compare(end1, end2);
      if (res != 0) {
        return res;
      }
      return super.compare(o1, o2);
    }
  }
}
