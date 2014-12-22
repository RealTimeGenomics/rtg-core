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

package com.rtg.variant.bayes.multisample;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import com.rtg.variant.Variant;
import com.rtg.variant.Variant.VariantFilter;
import com.rtg.variant.bayes.multisample.ComplexRegion.RegionType;

/**
 */
public final class OutputUtils {
  private OutputUtils() { }

  private static <T> T step(final Iterator<T> itA) {
    return itA.hasNext() ? itA.next() : null;
  }

  /**
   * @param calls original SNP (ie non-complex calls).
   * @param cx the complex regions extracted from calls.
   * @return the calls in <code>calls</code> that are not included in any of the non-hyper regions in <code>cx</code>.
   */
  static <T extends Variant> List<T> nonShadowed(List<T> calls, Complexities cx) {
    final Iterator<T> ita = calls.iterator();
    final List<T> res = new LinkedList<>();
    final Iterator<ComplexRegion> itb = cx.iterator();
    ComplexRegion cr = step(itb);
    while (ita.hasNext()) {
      final T call = ita.next();
      if (call.isIndel() || call.isOverflow()) {
        continue;
      }
      while (cr != null && cr.getEnd() <= call.getLocus().getStart()) {
        cr = step(itb);
      }
      if (cr == null) {
        res.add(call);
        continue;
      }
      if (call.getLocus().getStart() >= cr.getStart() && call.getLocus().getEnd() <= cr.getEnd()) {
        // System.err.println("call.getStart()=" + call.getLocus().getStart() + " call.getEnd()=" + call.getLocus().getEnd() + " cx.startPosition()=" + cr.startPosition() + " cr.endPosition()=" + cr.endPosition() + " " + cr.type());
        if (cr.type() == RegionType.HYPER) {
          //System.err.println("Failed region: " + cr);
          call.addFilter(VariantFilter.HYPER_COMPLEX);
          res.add(call);
        }
        continue;
      }
      assert cr.getStart() >= call.getLocus().getEnd();
      res.add(call);
    }
    return res;
  }

  /**
   * Merge two disjoint lists of variant calls.
   * @param callsA first list to be merged.
   * @param callsB second list to be merged.
   * @return merged list.
   */
  static <T extends Variant> List<T> merge(List<T> callsA, List<T> callsB) {
    final List<T> res = new LinkedList<>();
    final Iterator<T> itA = callsA.iterator();
    T callA = itA.hasNext() ? itA.next() : null;
    final Iterator<T> itB = callsB.iterator();
    T callB = itB.hasNext() ? itB.next() : null;
    while (callA != null || callB != null) {
      if (callA == null || (callB != null && callA.getLocus().getStart() > callB.getLocus().getStart())) {
        res.add(callB);
        callB = step(itB);
      } else if (callB == null || (callA.getLocus().getStart() < callB.getLocus().getStart())) { //we already know callA is non-null since it got past previous if
        res.add(callA);
        callA = step(itA);
      } else {
        assert callA.getLocus().getStart() == callB.getLocus().getStart();
        // Two calls have the same start position.  One of them should be length 0 for indel call
        if (callA.isOverflow()) {
          callA = step(itA);
        } else if (callB.isOverflow()) {
          callB = step(itB);
        } else if (callA.getLocus().getEnd() == callA.getLocus().getStart()) {
          res.add(callA);
          callA = step(itA);
        } else if (callB.getLocus().getEnd() == callB.getLocus().getStart()) {
          res.add(callB);
          callB = step(itB);
        } else {
          res.add(callA);
          callA = step(itA);

          // Two calls same start with non-trivial lengths
          //System.err.println(callA);
          //System.err.println(callB);
        }
      }
    }
    return res;
  }

}
