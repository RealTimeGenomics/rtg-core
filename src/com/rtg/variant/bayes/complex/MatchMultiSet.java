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
package com.rtg.variant.bayes.complex;

import static com.rtg.util.StringUtils.LS;

import java.util.Locale;
import java.util.Map;
import java.util.TreeMap;

import com.rtg.util.Utils;
import com.rtg.variant.match.Match;

/**
 * Holds a multiset with counts for each match.
 */
public class MatchMultiSet {

  private final Map<String, SingleCounts> mMatchMap = new TreeMap<>();
  private long mTotalCount = 0;

  /**
   * Construct an empty set.
   */
  public MatchMultiSet() { }

  MatchMultiSet(final String name) {
    mMatchMap.put(name, new SingleCounts());
  }

  /**
   * Add a new match.
   * @param match to be added.
   * @param corr correction value - needed as separate parameter because requires access to default q-score.
   */
  public void add(final Match match, final double corr) {
    final String hyp = match.toString();
    final SingleCounts count;
    if (mMatchMap.containsKey(hyp)) {
      count = mMatchMap.get(hyp);
    } else {
      count = new SingleCounts();
      mMatchMap.put(hyp, count);
    }
    count.increment(corr);
    ++mTotalCount;
  }

  /**
   * The number of distinct names for matches.
   * @return the number of names.
   */
  int size() {
    return mMatchMap.size();
  }

  /**
   * The total number of matches represented by this set.
   * @return total count of matches.
   */
  long totalCount() {
    return mTotalCount;
  }

  /**
   * Get all the distinct names of matches.
   * @return an array of names.
   */
  String[] names() {
    final int size = mMatchMap.size();
    return mMatchMap.keySet().toArray(new String[size]);
  }

  /**
   * @param lowerCase true if the output should be forced to lowercase
   * @return human readable output in the official output format
   */
  public String output(final boolean lowerCase) {
    final StringBuilder sb = new StringBuilder();
    for (final Map.Entry<String, SingleCounts> entry : mMatchMap.entrySet()) {
      final String name = entry.getKey();
      final SingleCounts count = entry.getValue();
      if (count.count() == 0) {
        continue;
      }
      sb.append('\t');
      final String na = lowerCase ? name.toLowerCase(Locale.getDefault()) : name;
      sb.append(na).append('\t').append(count.count()).append('\t').append(Utils.realFormat(count.correction(), 3));
    }
    return sb.toString();
  }

  /**
   * Reset the object
   */
  public void reset() {
    mMatchMap.clear();
    mTotalCount = 0;
  }

  @Override
  public String toString() {
    final StringBuilder sb = new StringBuilder();
    sb.append("[ MatchMultiSet size=").append(size()).append(LS);
    for (final Map.Entry<String, SingleCounts> entry : mMatchMap.entrySet()) {
      final String name = entry.getKey();
      final SingleCounts count = entry.getValue();
      sb.append(name).append(" > ").append(count.toString()).append(LS);
    }
    sb.append("]").append(LS);
    return sb.toString();
  }
}
