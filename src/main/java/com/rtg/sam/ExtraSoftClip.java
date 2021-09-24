/*
 * Copyright (c) 2016. Real Time Genomics Limited.
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

import java.util.List;

import com.rtg.util.diagnostic.Diagnostic;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

/**
 * Look for cases that were not soft-clipped during primary alignment but perhaps should have been.
 */
public final class ExtraSoftClip {

  /** If an alignment contains XID within this many bases of the end of the read, soft clip back to next match */
  private static final int MIN_ANCHOR = 3;

  private ExtraSoftClip() { }

  /**
   * Masks mismatches and inserts present at the end of the strand in given read
   * @param record sam record to modify
   * @return true if record was modified
   */
  public static boolean addSoftClip(SAMRecord record) {
    if (record.getReadUnmappedFlag()) {
      return false;
    }
    if (record.getReadNegativeStrandFlag()) {
      return checkNeg(record);
    } else {
      return checkPos(record);
    }
  }

  static boolean checkNeg(SAMRecord record) {
    int matchCount = 0;
    int scCount = 0;
    int mmCount = 0;
    final List<CigarElement> cigarElements = record.getCigar().getCigarElements();
    for (int l = 0; l < cigarElements.size(); l++) {
      final CigarElement e = cigarElements.get(l);
      final CigarOperator op = e.getOperator();
      switch (op) {
        case EQ:
          matchCount += e.getLength();
          if (matchCount >= MIN_ANCHOR) {
            if (mmCount == 0) {
              return false;
            }
            final Cigar c = new Cigar();
            if (scCount > 0) {
              c.add(new CigarElement(scCount, CigarOperator.S));
              record.setAlignmentStart(record.getAlignmentStart() + scCount);
            }
            do {
              c.add(cigarElements.get(l++));
            } while (l < cigarElements.size());
            record.setCigar(c);
            return true;
          }
          scCount += e.getLength();
          break;
        case X:
        case I:
        case D:
          if (op.consumesReadBases()) {
            scCount += e.getLength();
          }
          mmCount += e.getLength();
          break;
        default:
          if (l > 0 && matchCount <= MIN_ANCHOR) {
            Diagnostic.developerLog("Unhandled neg record " + matchCount + " -- " + record.getSAMString());
          }
          return false;
      }
    }
    if (matchCount <= MIN_ANCHOR) {
      //System.err.println("Setting record to unmapped with insufficient matches " + matchCount + " -- "  + record.getSAMString());
      SamUtils.convertToUnmapped(record);
      return true;
    }
    return false;
  }

  static boolean checkPos(SAMRecord record) {
    int matchCount = 0;
    int scCount = 0;
    int mmCount = 0;
    final List<CigarElement> cigarElements = record.getCigar().getCigarElements();
    for (int l = cigarElements.size() - 1; l >= 0; l--) {
      final CigarElement e = cigarElements.get(l);
      final CigarOperator op = e.getOperator();
      switch (op) {
        case EQ:
          matchCount += e.getLength();
          if (matchCount >= MIN_ANCHOR) {
            if (mmCount == 0) {
              return false;
            }
            final Cigar c = new Cigar();
            int i = 0;
            while (i <= l) {
              c.add(cigarElements.get(i++));
            }
            if (scCount > 0) {
              c.add(new CigarElement(scCount, CigarOperator.S));
            }
            record.setCigar(c);
            return true;
          }
          scCount += e.getLength();
          break;
        case X:
        case I:
        case D:
          if (op.consumesReadBases()) {
            scCount += e.getLength();
          }
          mmCount += e.getLength();
          break;
        default:
          if (l < cigarElements.size() - 1 && matchCount <= MIN_ANCHOR) {
            Diagnostic.developerLog("Unhandled pos record " + matchCount + " -- "  + record.getSAMString());
          }
          return false;
      }
    }
    if (matchCount <= MIN_ANCHOR) {
      //System.err.println("Setting record to unmapped with insufficient matches " + matchCount + " -- "  + record.getSAMString());
      SamUtils.convertToUnmapped(record);
      return true;
    }
    return false;
  }
}
