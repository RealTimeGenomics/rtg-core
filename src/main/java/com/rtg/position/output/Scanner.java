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
