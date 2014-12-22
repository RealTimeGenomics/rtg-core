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
package com.rtg.ngs.tempstage;

import java.io.Serializable;
import java.util.Comparator;

/**
 * Comparator of temp file records for Paired end data
 */
public class PairedTempFileRecordComparator implements Comparator<BinaryTempFileRecord>, Serializable {

  @Override
  public int compare(BinaryTempFileRecord o1, BinaryTempFileRecord o2) {
    assert o1.getReferenceId() == o2.getReferenceId();
    final int startCompare = o1.getStartPosition() - o2.getStartPosition();
    if (startCompare == 0) {
      final int readNameCompare = o1.getReadId() - o2.getReadId();
      if (readNameCompare != 0) {
        return readNameCompare < 0 ? -1 : 1;
      }
      final int mateStart = o1.getMatePosition() - o2.getMatePosition();
      if (mateStart == 0) {
        final boolean strandCompareDifferent = o1.isReverseStrand() ^ o2.isReverseStrand();
        if (!strandCompareDifferent) {
          final boolean pairedDifferent = o1.isReadPaired() ^ o2.isReadPaired();
          if (pairedDifferent) {
            return o1.isReadPaired() ? 1 : -1;
          } else {
            if (!o1.isReadPaired()) {
              return compareScores(o1, o2);
            }
          }
          final boolean mateStrandCompareDifferent = o1.isMateReverseStrand() ^ o2.isMateReverseStrand();
          if (!mateStrandCompareDifferent) {
            if (o1.isFirstOfPair() != o2.isFirstOfPair()) {
              return o1.isFirstOfPair() ? -1 : 1;
            }
            return compareScores(o1, o2);
          } else {
            return o1.isMateReverseStrand() ? -1 : 1;
          }
        }
        return o1.isReverseStrand() ? -1 : 1;
      } else {
        return mateStart;
      }

    } else {
      return startCompare;
    }
  }

  int compareScores(BinaryTempFileRecord o1, BinaryTempFileRecord o2) {
    final int scoreCompare = o1.getAlignmentScore() - o2.getAlignmentScore();
    if (scoreCompare != 0) {
      return scoreCompare;
    }
    return o1.getComboScore() - o2.getComboScore();
  }
}
