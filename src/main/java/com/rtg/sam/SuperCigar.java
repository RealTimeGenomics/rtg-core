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

import java.util.Arrays;

import com.rtg.alignment.ActionsHelper;
import com.rtg.alignment.ActionsHelper.CommandIterator;
import com.rtg.mode.DnaUtils;

/**
 * An extension of SAM cigars to handle CG overlaps.
 * These cigars can have all the same commands as SAM cigars,
 * plus a 'B' command that means move backwards along the template.
 * If a super cigar contains no B commands, then it is a normal SAM cigar.
 *
 * However, a super cigar usually goes hand-in-hand with a 'read-delta' string,
 * which contains just the nucleotides in the read that are different from the
 * template or missing from the template.  So a super cigar plus a read-delta
 * string are sufficient to reconstruct a read, given a template and a starting position.
 *
 */
public final class SuperCigar {

  /** Super CIGAR value for OVERLAP IN READ */
  public static final char OVERLAP_IN_READ = SamUtils.CIGAR_OVERLAP_IN_READ;
  /** Super CIGAR value for UNKNOWN IN TEMPLATE */
  public static final char UNKNOWN_TEMPLATE = SamUtils.CIGAR_UNKNOWN_TEMPLATE;
  /** Super CIGAR value for UNKNOWN IN READ - note if both template AND read are unknown, you should use this value. */
  public static final char UNKNOWN_READ = SamUtils.CIGAR_UNKNOWN_READ;

  private SuperCigar() {
  }

  /** Maps action numbers to super cigar characters (a superset of SAM cigars). */
  static final char[] SAM_CIGAR = new char[ActionsHelper.getNumActions()];
  static {
    Arrays.fill(SAM_CIGAR, (char) -1);
    SAM_CIGAR[ActionsHelper.SAME] = SamUtils.CIGAR_SAME;
    SAM_CIGAR[ActionsHelper.MISMATCH] = SamUtils.CIGAR_MISMATCH;
    SAM_CIGAR[ActionsHelper.DELETION_FROM_REFERENCE] = SamUtils.CIGAR_DELETION_FROM_REF;
    SAM_CIGAR[ActionsHelper.INSERTION_INTO_REFERENCE] = SamUtils.CIGAR_INSERTION_INTO_REF;
    SAM_CIGAR[ActionsHelper.CG_GAP_IN_READ] = SamUtils.CIGAR_GAP_IN_READ;
    SAM_CIGAR[ActionsHelper.CG_OVERLAP_IN_READ] = OVERLAP_IN_READ;
    SAM_CIGAR[ActionsHelper.SOFT_CLIP] = SamUtils.CIGAR_SOFT_CLIP;
    SAM_CIGAR[ActionsHelper.UNKNOWN_READ] = UNKNOWN_READ;
    SAM_CIGAR[ActionsHelper.UNKNOWN_TEMPLATE] = UNKNOWN_TEMPLATE;
  }

  private static void updateCigar(StringBuilder sb, int action, int count, boolean softClip) {
    if (count > 0) {
        sb.append(count);
        sb.append(softClip && action != ActionsHelper.CG_OVERLAP_IN_READ ? SamUtils.CIGAR_SOFT_CLIP : SAM_CIGAR[action]);
    }
  }

  /**
   * Converts an actions array into a super cigar string.
   * This ignores clipping.
   *
   * @param actions packed actions array.
   * @param reverse the resulting cigar will be reversed if this is true.
   * @param templateLength total length of the template
   * @return super cigar
   */
  public static String actionsToSuperCigar(int[] actions, boolean reverse, int templateLength) {
    final CommandIterator iter = reverse ? ActionsHelper.iteratorReverse(actions) : ActionsHelper.iterator(actions);
    final StringBuilder sb = new StringBuilder();
    int prevAction = -1;
    int count = -1;
    long pos = actions[ActionsHelper.TEMPLATE_START_INDEX];
    while (iter.hasNext()) {
      final int action = iter.next();
      if (action != prevAction || (action != ActionsHelper.CG_OVERLAP_IN_READ && (pos == 0 || pos == templateLength))) {
        updateCigar(sb, prevAction, count, pos <= 0 || pos > templateLength);

        prevAction = action;
        count = 1;
      } else {
        ++count;
      }
      if (prevAction != ActionsHelper.CG_OVERLAP_IN_READ) {
        ++pos;
      } else {
        --pos;
      }
    }
    updateCigar(sb, prevAction, count, pos <= 0 || pos > templateLength);
    return sb.toString();
  }

  /**
   * Returns the read nucleotides that are different from the template.
   * This string, together with the super cigar, allow the original read
   * to be reconstructed from the template.
   *
   * @param actions packed actions array.
   * @param read the original read sequence.
   * @param templateLength total length of the template
   * @return a string containing just the missing and modified read nucleotides.
   */
  public static String readDelta(int[] actions, byte[] read, int templateLength) {
    final CommandIterator iter = ActionsHelper.iterator(actions);
    final StringBuilder sb = new StringBuilder();
    int readPos = 0;
    int pos = actions[ActionsHelper.TEMPLATE_START_INDEX];

    while (iter.hasNext()) {
      final int action = iter.next();
      //System.err.println(pos + " : " + templateLength + " : " + readPos + " :: " + SAM_CIGAR[action]);
      if (action == ActionsHelper.MISMATCH || action == ActionsHelper.INSERTION_INTO_REFERENCE || action == ActionsHelper.UNKNOWN_TEMPLATE
          || (action != ActionsHelper.CG_OVERLAP_IN_READ && action != ActionsHelper.CG_GAP_IN_READ && (pos < 0 || pos >= templateLength))) {
        sb.append(DnaUtils.getBase(read[readPos]));
        ++readPos;
        ++pos;
      } else if (action == ActionsHelper.SAME || action == ActionsHelper.UNKNOWN_READ) {  //read is known
        ++readPos;
        ++pos;
      } else  if (action == ActionsHelper.CG_GAP_IN_READ) {
        ++pos;
      } else if (action == ActionsHelper.CG_OVERLAP_IN_READ) {
        --pos;
      }
    }
    assert readPos == read.length;
    return sb.toString();
  }
}
