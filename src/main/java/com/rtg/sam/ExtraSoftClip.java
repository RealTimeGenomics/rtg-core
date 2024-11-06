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
