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
  static List<Variant> nonShadowed(List<Variant> calls, Complexities cx) {
    final Iterator<Variant> ita = calls.iterator();
    final List<Variant> res = new LinkedList<>();
    final Iterator<ComplexRegion> itb = cx.iterator();
    ComplexRegion cr = step(itb);
    while (ita.hasNext()) {
      final Variant call = ita.next();
      if (call.isIndel() || call.isOverflow() || call.isSoftClip()) {
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
  static List<Variant> merge(List<Variant> callsA, List<Variant> callsB) {
    final List<Variant> res = new LinkedList<>();
    final Iterator<Variant> itA = callsA.iterator();
    Variant callA = itA.hasNext() ? itA.next() : null;
    final Iterator<Variant> itB = callsB.iterator();
    Variant callB = itB.hasNext() ? itB.next() : null;
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
