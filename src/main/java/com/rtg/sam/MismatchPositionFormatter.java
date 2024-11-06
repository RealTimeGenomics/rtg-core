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
      for (int i = zeroBasedStart; i < 0 && iter.hasNext(); ++i) {
        iter.next();
      }
    }
    int templatePosition = Math.max(0, zeroBasedStart);
    boolean inDeletion = false;
    int matchCount = 0;
    while (templatePosition < template.length && iter.hasNext()) {
      final int action = iter.next();
      switch (action) {
        case ActionsHelper.SAME:
          ++matchCount;
          inDeletion = false;
          ++templatePosition;
          break;
        case ActionsHelper.MISMATCH:
          mismatch.append(matchCount);
          matchCount = 0;
          mismatch.append(DnaUtils.getBase(template[templatePosition]));
          inDeletion = false;
          ++templatePosition;
          break;
        case ActionsHelper.DELETION_FROM_REFERENCE:
          if (!inDeletion) {
            mismatch.append(matchCount);
            mismatch.append("^");
            matchCount = 0;
          }
          mismatch.append(DnaUtils.getBase(template[templatePosition]));
          inDeletion = true;
          ++templatePosition;
          break;
        case ActionsHelper.CG_GAP_IN_READ: // Probably not set up correctly
          ++templatePosition;
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
