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

package com.rtg.segregation;

import java.util.Deque;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map;
import java.util.NavigableSet;
import java.util.TreeSet;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.Pair;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Search for the lowest cost path through a sequence of blocks.
 */
//TODO some serious testing especially of garbage collection
@TestClass("com.rtg.segregation.SegregationVcfSearchTest")
public class Search extends IntegralAbstract {

  private static final int ERROR_PENALTY = 1;

  private static final int REMOVE_THRESHOLD = 10000;

  private final int mNewPenalty;

  private final int mXOPenalty;

  static final int DEFAULT_NEW_PENALTY = 1000;
  static final int DEFAULT_XO_PENALTY = 100;

  Search(int newPenalty, int xOPenalty) {
    mNewPenalty = newPenalty;
    mXOPenalty = xOPenalty;
  }

  private int mSearchContainderId = 0;

  private  NavigableSet<SearchContainer> mRankedChains = new TreeSet<>();

  void add(final SegregationBlock block) {
    //System.err.println("add " + block.toString());
    final Map<Pair<PatternArray, Boolean>, SearchContainer> map = new HashMap<>();
    //Add a new break (possibly the first one) - use the best chain so far.
    final SearchContainer n0;
    final PatternArray patterns = block.patterns();
    if (mRankedChains.isEmpty()) {
      final SearchContainer n1 = new SearchContainer(null, block.count() * ERROR_PENALTY, null, block, patterns, SearchType.Error, null, mSearchContainderId);
      add(map, n1);
      n0 = new SearchContainer(null, 0, null, block, patterns, SearchType.New, null, mSearchContainderId);
    } else {
      final SearchContainer bestSoFar = mRankedChains.first();
      n0 = new SearchContainer(bestSoFar, bestSoFar.score() + mNewPenalty, null, block, patterns, SearchType.New, null, mSearchContainderId);
    }
    add(map, n0);

    //Check everything so far.
    for (final SearchContainer sc : mRankedChains) {
      //System.err.println(">>" + sc.toString());
      final SearchContainer good = sc.goodContainer();
      if (good == null) {
        final SearchContainer n1 = new SearchContainer(sc, sc.score(), null, block, patterns, SearchType.OK, null, mSearchContainderId);
        add(map, n1);
        continue;
      }
      //If last ok is cross over then do that
      final CrossOver crossover;
      if (sc.goodContainer().block().isXLike() != block.isXLike()) {
        crossover = null;
      } else {
        crossover = PatternArray.crossover(sc.pattern(), patterns, block.isXLike());
      }
      if (crossover != null) {
        //System.err.println("before=" + sc.pattern());
        //System.err.println("after=" + patterns);
        //System.err.println("cross=" + crossover);
        //final SearchContainer n1 = new SearchContainer(sc, sc.score() + mXOPenalty, null, block,  patterns, SearchType.XO, crossover, mSearchContainderId);
        final SearchContainer n1 = new SearchContainer(sc, sc.score() + mXOPenalty, null, block, crossover.pattern(), SearchType.XO, crossover, mSearchContainderId);
        add(map, n1);
      }

      //If good is compatible then add an ok for each possible flip
      //else add an error

      boolean ok = false;
      final PatternArray pa = good.pattern();
      for (int flip = 0; flip < Pattern.NUMBER_FLIPS; ++flip) {
        final PatternArray fi = pa.flipIntersect(patterns, flip);
        if (fi == null) {
          continue;
        }
        final SearchContainer n1 = new SearchContainer(sc, sc.score(), null, block, fi, SearchType.OK, null, mSearchContainderId);
        add(map, n1);
        ok = true;
      }

      if (!ok) {
        //add an error case.
        final SearchContainer n1 = new SearchContainer(sc, sc.score() + block.count() * ERROR_PENALTY, good, block, pa, SearchType.Error, null, mSearchContainderId);
        add(map, n1);
      }
    }

    //garbage collect entries when too many or too high a score
    final NavigableSet<SearchContainer> newChains = new TreeSet<>(map.values());
    final double removeScore = newChains.first().score() + REMOVE_THRESHOLD;
    final Iterator<SearchContainer> it = newChains.descendingIterator();
    int remaining = newChains.size();
    while (it.hasNext()) {
      final SearchContainer next = it.next();
      if (remaining > 1000 || next.score() > removeScore) {
        //System.err.println("##" + next.toString());
        it.remove();
        --remaining;
      } else {
        break;
      }
    }

    //Handy to turn on when debugging the search
    //    final Iterator<SearchContainer> it0 = newChains.iterator();
    //    while (it0.hasNext()) {
    //      final SearchContainer next = it0.next();
    //      //System.err.println("&&" + next.toString());
    //    }

    mRankedChains = newChains;
  }

  private void add(final Map<Pair<PatternArray, Boolean>, SearchContainer> map, final SearchContainer n0) {
    //System.err.println("++ " + n0.toString());
    final Pair<PatternArray, Boolean> key = n0.key();
    final SearchContainer earlier = map.get(key);
    final SearchContainer ni;
    if (earlier != null) {
      if (earlier.score() <= n0.score()) {
        //best is already there do nothing
        //System.err.println("-- " + n0.toString());
        ni = null;
      } else {
        //System.err.println("-- " + earlier.toString());
        ni = n0;
      }
    } else {
      ni = n0;
    }
    if (ni != null) {
      map.put(n0.key(), n0);
    }
    ++mSearchContainderId;
  }

  /**
   * @return iterator over the search containers for the best path in forward order.
   */
  Iterable<SearchContainer> bestResult() {
    if (mRankedChains.isEmpty()) {
      return null;
    }
    final SearchContainer best = mRankedChains.first();
    //scan chain in reverse
    final Deque<SearchContainer> que = new LinkedList<>();
    SearchContainer x = best;
    while (true) {
      que.addFirst(x);
      final SearchContainer next = x.previous();
      if (next == null) {
        break;
      }
      x = next;
    }
    return que;
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mSearchContainderId >= 0);
    Exam.assertNotNull(mRankedChains);
    return true;
  }

}
