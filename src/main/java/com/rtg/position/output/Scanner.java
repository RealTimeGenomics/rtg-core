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
package com.rtg.position.output;

/**
 * Placed in a separate class so the logic of scanning over buckets
 * can be tested.
 */
abstract class Scanner <G> {

  /**
   * Do scanning over all necessary buckets.
   * @param region passed through to underlying scan.
   * @param forward flag that says whether looking for matching regions forward from region.
   * @param first bucket associated with region.
   * @param last the last bucket associated with the region (inclusive).
   * @param lo start position for scan.
   * @param hi the last position for the scan (inclusive).
   * @return the best matching region.
   */
  G scanAll(final G region, final boolean forward, final long first, final long last, final long lo, final long hi) {
    init();
    //System.err.println("scanall first=" + first + " last=" + last + " lo=" + lo + " hi=" + hi);
    assert first <= last && lo <= hi;
    final long len = last - first + 1;
    final long les = hi - lo + 1;
    if (les >= len) {
      //scan everything
      scan(region, forward, first, last);
    } else if (hi < first) {
      final long l = first - hi;
      final long m = (l - 1) / len + 1;
      assert m >= 1;
      final long ml = m * len;
      final long lo1 = lo + ml;
      final long hi1 = hi + ml;
      assert first <= hi1;
      assert lo1 <= last;
      assert lo1 <= hi1;
      scanAll2(region, forward, first, last, lo1, hi1);
    } else if (lo > last) {
      final long l = lo - last;
      final long m = (l - 1) / len + 1;
      final long ml = m * len;
      final long lo1 = lo - ml;
      assert lo1 <= last;
      final long hi1 = hi - ml;
      assert hi1 >= first;
      assert lo1 <= hi1;
      scanAll2(region, forward, first, last, lo1, hi1);
    } else {
      scanAll2(region, forward, first, last, lo, hi);
    }
    return result();
  }


  private void scanAll2(final G region, final boolean forward, final long first, final long last, final long lo, final long hi) {
    assert lo <= hi; // : lo + "=lo " + hi + "=hi";
    if (first <= lo && hi <= last) {
      scan(region, forward, lo, hi);
    } else if (hi > last) {
      assert lo > first && lo <= last;
      final long hi1 = hi - last + first - 1;
      assert hi1 < lo;
      assert first <= hi1;
      scan(region, forward, first, hi1);
      assert lo <= last;
      if (resultScore() != 0.0) {
        scan(region, forward, lo, last);
        assert (hi - lo + 1) == (hi1 - first + 1) + (last - lo + 1);
      }
    } else if (lo < first) {
      assert hi < last && first <= hi;
      final long lo1 = lo - first + last + 1;
      assert lo1 > hi;
      assert first <= hi;
      scan(region, forward, first, hi);
      assert lo1 <= last;
      if (resultScore() != 0.0) {
        scan(region, forward, lo1, last);
        assert (hi - lo + 1) == (hi - first + 1) + (last - lo1 + 1);
      }
    } else {
      throw new IllegalArgumentException();
    }
  }

/**
   * Scan the buckets from lo to hi (inclusive).
   * @param region to be checked when scanning.
   * @param forward true if doing a forward scan.
   * @param lo first bucket to be scanned.
   * @param hi last bucket to be scanned (inclusive).
   */
  protected abstract void scan(final G region, final boolean forward, final long lo, final long hi);

  /** Called once before start of 1 or 2 scans. */
  protected abstract void init();

  /**
   * Get the best match region found by the scans.
   * @return the best region from the underlying scans (may be null if no matches found).
   */
  protected abstract G result();

  /**
   * Return the score for the best match region found by the scans.
   * @return the score for the best region from the underlying scans.

   */
  protected abstract double resultScore();
}
