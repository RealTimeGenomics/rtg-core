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

package com.rtg.variant.bayes.complex;

import com.rtg.reference.Ploidy;
import com.rtg.util.StringUtils;
import com.rtg.variant.Variant;
import com.rtg.variant.util.VariantUtils;

/**
 * Filter out calls exhibiting probably homopolymer problems.
 */
public final class IonTorrentCallFilter {

  private IonTorrentCallFilter() { }

  /**
   * Determine whether this call should be filtered due to detected Ion Torrent homopolymer badness
   * @param template the template
   * @param res the variant call to check
   * @param useFilter whether to use this filter or not
   * @return true if the call should be filtered
   */
  public static boolean ionTorrentFilter(byte[] template, Variant res, boolean useFilter) {
    if (!useFilter) {
      return false;
    }

    final String ref = res.getLocus().getRefNts();
    final String name = res.getSample(0).getName();
    if (name == null) {
      return false;
    }

    final String call;
    if (res.getSample(0).getPloidy() == Ploidy.DIPLOID) {  //TODO polyploid will suck!
      final String[] calls = StringUtils.split(name, VariantUtils.COLON);
      if (calls[0].equals(calls[1])) {
        call = calls[0];
      } else {
        final boolean firstEqRef = calls[0].equals(ref);
        final boolean secondEqRef = calls[1].equals(ref);
        if (!firstEqRef && !secondEqRef) {
          return false; //heterozygous and neither side are the reference, just assume it's not homopolymer
        }
        call = firstEqRef ? calls[1] : calls[0];
      }
    } else {
      call = name;
    }

    if (ref.length() == call.length()) { // not an indel...
      return false;
    }
    final char ch0 = ref.length() > 0 ? ref.charAt(0) : call.charAt(0);
    for (int i = 1; i < ref.length(); ++i) {
      if (ref.charAt(i) != ch0) {
        return false;
      }
    }
    for (int i = 0; i < call.length(); ++i) {
      if (call.charAt(i) != ch0) {
        return false;
      }
    }
    int incrementD = 0;
    final int start = res.getLocus().getStart();
    final int nt0 = template[start];
    int i = start - 1;
    while (i >= 0 && nt0 == template[i]) {
      --i;
      ++incrementD;
    }
    int incrementU = 0;
    int j = start + ref.length();
    while (j < template.length && nt0 == template[j]) {
      ++j;
      ++incrementU;
    }
    final int increment = incrementD + incrementU;
    return ionFilter(nt0, ref.length() + increment, increment + call.length());
  }

  private static final boolean[][][] IONT_FILTER = new boolean[5][5][5];
  static {
    IONT_FILTER[1][2][1] = true;  //AA->A
    IONT_FILTER[2][2][1] = true;  //CC->C
    IONT_FILTER[3][2][1] = true;  //GG->G
    IONT_FILTER[2][3][2] = true;  //CCC->CC
    IONT_FILTER[3][3][2] = true;  //GGG->GG
    IONT_FILTER[4][3][4] = true;  //TTT->TTTT
    IONT_FILTER[2][4][3] = true;  //CCCC->CCC
  }

  static boolean ionFilter(int nt0, int i, int j) {
    if (nt0 >= IONT_FILTER.length) {
      return false;
    }
    final boolean[][] ix1 = IONT_FILTER[nt0];
    if (i >= ix1.length) {
      return false;
    }
    final boolean[] ix2 = ix1[i];
    if (j >= ix2.length) {
      return false;
    }
    return ix2[j];
  }
}
