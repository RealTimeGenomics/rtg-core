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

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

/**
 * Detect and mask homopolymer mismatches near ends of read
 */
public final class HomopolymerUtils {

  private static final int MIN_HOMO = 3;

  private HomopolymerUtils() { }

  /**
   * Masks mismatches and inserts present at the end of the strand in given read
   * @param record sam record to modify
   * @return true if record was modified
   */
  public static boolean maskHomoPolymer(SAMRecord record) {
    if (record.getReadNegativeStrandFlag()) {
      return checkNeg(record);
    } else {
      return checkPos(record);
    }
  }

  static boolean checkNeg(SAMRecord record) {
    boolean modified = false;
    int readPos = 0;
    final byte[] readBases = record.getReadBases();
    boolean isHomo = true;
    byte homoBase = -1;
    if (readBases.length > 0) {
      homoBase = readBases[0];
    }
    for (int i = 0; i < MIN_HOMO && i < readBases.length; ++i) {
      if (readBases[i] != homoBase) {
        isHomo = false;
      }
    }
    for (CigarElement e : record.getCigar()) {
      final CigarOperator op = e.getOperator();
      final int length = e.getLength();
      if (isHomo) {
        if (op.consumesReadBases()) {
          for (int i = readPos; i < readPos + length; ++i) {
            if (readBases[i] == homoBase) {
              modified = maskBase(readBases, op, i) || modified;
            } else {
              isHomo = false;
              break;
            }
          }
          readPos += length;
        }
      } else {
        break;
      }
    }
    record.setReadBases(readBases);
    return modified;
  }

  static boolean checkPos(SAMRecord record) {
    boolean modified = false;
    final byte[] readBases = record.getReadBases();
    int readPos = readBases.length - 1;
    boolean isHomo = true;
    byte homoBase = -1;
    if (readBases.length > 0) {
      homoBase = readBases[readPos];
    }
    for (int i = readPos; i > readPos - MIN_HOMO && i >= 0; --i) {
      if (readBases[i] != homoBase) {
        isHomo = false;
      }
    }

    final List<CigarElement> cigarElements = record.getCigar().getCigarElements();
    for (int l = cigarElements.size() - 1; l >= 0; --l) {
      final CigarElement e = cigarElements.get(l);
      final CigarOperator op = e.getOperator();
      final int length = e.getLength();
      if (isHomo) {
        if (op.consumesReadBases()) {
          for (int i = readPos; i > readPos - length && i >= 0; --i) {
            if (readBases[i] == homoBase) {
              modified = maskBase(readBases, op, i) || modified;
            } else {
              isHomo = false;
              break;
            }
          }
          readPos -= length;
        }
      } else {
        break;
      }
    }
    record.setReadBases(readBases);
    return modified;
  }

  private static boolean maskBase(byte[] readBases, CigarOperator op, int i) {
    if (op == CigarOperator.X  || op == CigarOperator.I) {
      readBases[i] = 'N';
      return true;
    }
    return false;
  }
}
