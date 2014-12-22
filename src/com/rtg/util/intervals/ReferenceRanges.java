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

package com.rtg.util.intervals;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;

import com.rtg.util.diagnostic.NoTalkbackSlimException;

/**
 * Contains a mapping from reference sequence name to RangeList. Ideally this would be merged with ReferenceRegions.
 */
public class ReferenceRanges {

  private final boolean mAllAvailable;

  private final SortedMap<String, RangeList<String>> mByName = new TreeMap<>();
  private SortedMap<Integer, RangeList<String>> mById = null;

  /**
   * Constructor
   * @param allAvailable true if satisfying these ranges should be done by streaming all data through (i.e. unrestricted).
   */
  public ReferenceRanges(boolean allAvailable) {
    mAllAvailable = allAvailable;
  }

  // Convert a map keyed off sequence names into one keyed off sequence ids using the supplied lookup
  private static <T> SortedMap<Integer, T> convertNameToIdKeys(Map<String, Integer> lookup, Map<String, T> strRanges) {
    // Convert from sequence name keys to sequenceId
    final SortedMap<Integer, T> rangeMap = new TreeMap<>();
    for (Map.Entry<String, T> entry : strRanges.entrySet()) {
      final Integer sequenceId = lookup.get(entry.getKey());
      if (sequenceId == null) {
        throw new NoTalkbackSlimException("Sequence \"" + entry.getKey() + "\" referenced in regions not found in the sequence dictionary.");
      }
      rangeMap.put(sequenceId, entry.getValue());
    }
    return rangeMap;
  }


  /**
   * @return true if no restrictions are in place that would require indexed input.
   */
  public boolean allAvailable() {
    return mAllAvailable;
  }


  /**
   * Add the supplied range list associated with the given sequence
   * @param seqName the name of the reference sequence
   * @param range the range data for the sequence
   */
  public void put(String seqName, RangeList<String> range) {
    mByName.put(seqName, range);
    if (mById != null) {
      throw new IllegalStateException("Cannot call put after setting sequence ids");
    }
  }

  /**
   * Use the supplied sequence name to id mapping to create an id based lookup
   * @param sequenceLookup the mapping from names to ids
   */
  public void setIdMap(Map<String, Integer> sequenceLookup) {
    mById = convertNameToIdKeys(sequenceLookup, mByName);
  }

  // Accessing via sequence names

  /**
   * Gets the regions associated with a sequence name
   * @param seqName the reference sequence
   * @return the regions, or null if none
   */
  public RangeList<String> get(String seqName) {
    return mByName.get(seqName);
  }

  /**
   * @param seqName the reference sequence
   * @return true if there are regions associated with a sequence
   */
  public boolean containsSequence(String seqName) {
    return mByName.containsKey(seqName);
  }

  /**
   * @return the set of known reference sequence names
   */
  public Collection<String> sequenceNames() {
    return mByName.keySet();
  }


  // Accessing via sequence ids

  /**
   * Gets the regions associated with a sequence id
   * @param seqId the reference sequence
   * @return the regions, or null if none
   */
  public RangeList<String> get(Integer seqId) {
    return mById.get(seqId);
  }

  /**
   * @param seqId the reference sequence
   * @return true if there are regions associated with a sequence
   */
  public boolean containsSequence(Integer seqId) {
    return mById.containsKey(seqId);
  }

  /**
   * @return the set of known reference sequence ids
   */
  public Collection<Integer> sequenceIds() {
    return mById.keySet();
  }

  @Override
  public String toString() {
    return mByName.toString();
  }

  /**
   * A helper class that allows incrementally adding RangeData and then conversion to a ReferenceRanges.
   */
  public static class Accumulator extends TreeMap<String, List<RangeList.RangeData<String>>> {

    /**
     * Adds a range data element
     * @param sequenceName the sequence that the range element is associated with
     * @param rangeData the range data element
     */
    public void addRangeData(String sequenceName, RangeList.RangeData<String> rangeData) {
      if (rangeData != null) {
        final List<RangeList.RangeData<String>> annos;
        if (containsKey(sequenceName)) {
          annos = get(sequenceName);
        } else {
          annos = new ArrayList<>();
          put(sequenceName, annos);
        }
        annos.add(rangeData);
      }
    }

    /**
     * Merges overlaps between regions and returns the result as a map of <code>RangeList</code> objects, keyed by chromosome name.
     *
     * @return a map of <code>RangeList</code> objects keyed by chromosome name.
     */
    public ReferenceRanges getReferenceRanges() {
      final ReferenceRanges rangeLists = new ReferenceRanges(false);
      for (final Map.Entry<String, List<RangeList.RangeData<String>>> me : entrySet()) {
        final RangeList<String> search = new RangeList<>(me.getValue());
        rangeLists.put(me.getKey(), search);
      }
      return rangeLists;
    }
  }
}
