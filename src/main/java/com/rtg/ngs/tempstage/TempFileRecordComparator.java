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
package com.rtg.ngs.tempstage;

import java.io.Serializable;
import java.util.Comparator;

/**
 * Comparator to use for removing duplicate records.
 * It is assumed both records are on the same reference sequence.
 */
public class TempFileRecordComparator implements Comparator<BinaryTempFileRecord>, Serializable {

  @Override
  public int compare(BinaryTempFileRecord o1, BinaryTempFileRecord o2) {
    assert o1.getReferenceId() == o2.getReferenceId();
    final int startCompare = o1.getStartPosition() - o2.getStartPosition();
    if (startCompare == 0) {
      final int readNameCompare = o1.getReadId() - o2.getReadId();
      if (readNameCompare != 0) {
        return readNameCompare < 0 ? -1 : 1;
      }
      final boolean strandCompareDifferent = o1.isReverseStrand() ^ o2.isReverseStrand();
      if (strandCompareDifferent) {
        return o1.isReverseStrand() ? -1 : 1;
      }
      return compareScores(o1, o2);
    } else {
      return startCompare;
    }
  }

  int compareScores(BinaryTempFileRecord o1, BinaryTempFileRecord o2) {
    final int scoreCompare = o1.getAlignmentScore() - o2.getAlignmentScore();
    if (scoreCompare != 0) {
      return scoreCompare;
    }
    return 0;
  }
}
