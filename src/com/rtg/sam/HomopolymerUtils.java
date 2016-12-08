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
