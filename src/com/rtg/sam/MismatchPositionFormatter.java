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
package com.rtg.sam;

import com.rtg.alignment.ActionsHelper;
import com.rtg.mode.DnaUtils;

/**
 * This class is used to create a mismatching positions (MD) attribute string.
 */
public final class MismatchPositionFormatter {

  private MismatchPositionFormatter() { }

  /**
   * Converts an array of actions and a template array into a mismatching positions (MD) string.
   * Will not work for CG reads.
   * @param actions the alignment actions
   * @param rc is the read reverse complement
   * @param template template array
   * @return a mismatching positions (MD) string
   */
  public static String actionsToMismatchPositions(int[] actions, boolean rc, byte[] template) {
    final StringBuilder mismatch = new StringBuilder();
    final int zeroBasedStart = ActionsHelper.zeroBasedTemplateStart(actions);
    final ActionsHelper.CommandIterator iter = rc ? ActionsHelper.iteratorReverse(actions) : ActionsHelper.iterator(actions);
    if (zeroBasedStart < 0) {
      for (int i = zeroBasedStart; i < 0 && iter.hasNext(); i++) {
        iter.next();
      }
    }
    int templatePosition = Math.max(0, zeroBasedStart);
    boolean inDeletion = false;
    int matchCount = 0;
    while (iter.hasNext() && templatePosition < template.length) {
      final int action = iter.next();
      switch (action) {
        case ActionsHelper.SAME:
          matchCount++;
          inDeletion = false;
          templatePosition++;
          break;
        case ActionsHelper.MISMATCH:
          mismatch.append(matchCount);
          matchCount = 0;
          mismatch.append(DnaUtils.getBase(template[templatePosition]));
          inDeletion = false;
          templatePosition++;
          break;
        case ActionsHelper.DELETION_FROM_REFERENCE:
          if (!inDeletion) {
            mismatch.append(matchCount);
            mismatch.append("^");
            matchCount = 0;
          }
          mismatch.append(DnaUtils.getBase(template[templatePosition]));
          inDeletion = true;
          templatePosition++;
          break;
        case ActionsHelper.CG_GAP_IN_READ: // Probably not set up correctly
          templatePosition++;
          break;
        case ActionsHelper.CG_OVERLAP_IN_READ:
        case ActionsHelper.INSERTION_INTO_REFERENCE: // Ignore insertions on reference
        default:
          break;
      }
    }
    mismatch.append(matchCount);

    return mismatch.toString();
  }
}
