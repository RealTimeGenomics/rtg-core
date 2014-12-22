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

import com.rtg.reference.Ploidy;
import com.rtg.util.StringUtils;
import com.rtg.variant.Variant;

/**
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
    //    res.getName();
    //    res.getRef();
    //    res.getUntrimmedRef();
    //    res.getRef();
    //    System.err.print("name=" + res.getName() + " ref=" + res.getRef() + " untrimmedRef=" + res.getUntrimmedRef());
    //    System.err.print(" position ref:" + template[res.getPosition()] + "start ref:" + template[res.getStart()] + "output ref:" + template[res.getOutputPosition()]);
    //    for (int i = res.getStart(); i < res.getEnd(); i++) {
    //      System.err.print(template[res.getPosition()]);
    //    }
    //    System.err.println();
    //    res.getBestCat().name();

    final String ref = res.getLocus().getRefNts();
    final String name = res.getSample(0).getName();
    if (name == null) {
      return false;
    }

    final String call;
    if (res.getSample(0).getPloidy() == Ploidy.DIPLOID) {  //TODO polyploid will suck!
      final String[] calls = StringUtils.split(name, ':');
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
    for (int i = 1; i < ref.length(); i++) {
      if (ref.charAt(i) != ch0) {
        return false;
      }
    }
    for (int i = 0; i < call.length(); i++) {
      if (call.charAt(i) != ch0) {
        return false;
      }
    }
    int incrementD = 0;
    final int start = res.getLocus().getStart();
    final int nt0 = template[start];
    int i = start - 1;
    while (i >= 0 && nt0 == template[i]) {
      i--;
      incrementD++;
    }
    int incrementU = 0;
    int j = start + ref.length();
    while (j < template.length && nt0 == template[j]) {
      j++;
      incrementU++;
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
