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
package com.rtg.visualization;

import java.util.HashMap;
import java.util.Locale;

import com.rtg.sam.SamUtils;
import com.rtg.sam.SuperCigar;

import htsjdk.samtools.SAMRecord;

/**
 * This class is used for unrolling <code>CIGARS</code>
 *
 */
final class CigarHelper {

  private static final String INVALID_CIGAR = "Invalid cigar : ";

  private CigarHelper() {
  }

  static boolean isMissing(final String cigar) {
    return cigar == null || "*".equals(cigar);
  }

  /**
   * Given a cigar, locate inserts.
   *
   * @param cigar cigar to process
   * @param offset start
   * @param inserts insert locations
   */
  static void locateInserts(final String cigar, final int offset, final int[] inserts) {
    int n = 0;
    int rPos = 0;
    for (int i = 0; i < cigar.length(); ++i) {
      final char c = cigar.charAt(i);
      if (Character.isDigit(c)) {
        n = 10 * n + c - '0';
      } else {
        //System.err.println(n + " of " + c);
        assert n > 0;
        switch (c) {
          case SamUtils.CIGAR_INSERTION_INTO_REF:
            final int newOffset = rPos + offset;
            inserts[newOffset] = Math.max(inserts[newOffset], n);
            break;
          case SamUtils.CIGAR_DELETION_FROM_REF:
          case SamUtils.CIGAR_GAP_IN_READ:
          case SamUtils.CIGAR_SAME_OR_MISMATCH:
          case SamUtils.CIGAR_SAME:
          case SamUtils.CIGAR_MISMATCH:
          case SuperCigar.UNKNOWN_TEMPLATE:
          case SuperCigar.UNKNOWN_READ:
            rPos += n;
            break;
          case SamUtils.CIGAR_SOFT_CLIP:
          case SamUtils.CIGAR_HARD_CLIP:
            break;
          case 'B':
            rPos -= n;
            break;
          default:
            throw new IllegalStateException(INVALID_CIGAR + cigar);
        }
        n = 0;
      }
    }
  }


  static HashMap<Integer, String> unrollLegacyCgCigar(SAMRecord rec) {
    final HashMap<Integer, String> map = new HashMap<>();
    final String res = unrollCigar(rec);
    // System.err.println("unroll=\n" + res);

    final String gc = rec.getStringAttribute(SamUtils.ATTRIBUTE_CG_RAW_READ_INSTRUCTIONS);
    final String gs = rec.getStringAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_BASES);

    final StringBuilder sb = new StringBuilder();
    int n = 0;
    int readPos = 0;
    Integer startPos = 0;
    int overlapSize;
    int totalProcessed = 0;
    for (int i = 0; i < gc.length(); ++i) {
      final char c = gc.charAt(i);
      if (Character.isDigit(c)) {
        n = 10 * n + c - '0';
      } else {
        for (int j = 0; j < n; ++j) {
          ++totalProcessed;
          switch (c) {
            case 'G' :
              overlapSize = n;
              // System.err.println("overSize=" + overlapSize + " n=" + n + " sb.length()=" + sb.length());
              for (int k = overlapSize, m = 0; k > 0 ; --k, ++m) {
                if (totalProcessed > (rec.getReadLength() / 2)) {
                  sb.append(res.charAt(readPos++));
                } else {
                  sb.append(gs.charAt(m));
                }
              }
              //   System.err.println("kept in map = " + sb.toString() + "at pos = " + startPos);
              map.put(startPos, sb.toString());
              sb.setLength(0);
              for (int k = overlapSize; k > 0 ; --k) {
                if (totalProcessed > (rec.getReadLength() / 2)) {
                  sb.append(gs.charAt(overlapSize - k));
                } else {
                  sb.append(res.charAt(readPos++));
                }
              }
              startPos = readPos - overlapSize;
              j = n;
              break;
            case 'S' :
              //    System.err.println("j = " + j + " readPos=" + readPos);
              if (res.charAt(readPos) == ' ') {
                --j;
              }
              sb.append(res.charAt(readPos));
              ++readPos;
              break;
            default :
              throw new IllegalArgumentException("unknown code " + c);
          }
        }
        n = 0;
      }
    }
    map.put(startPos, sb.toString());
    //printmap(map);
    return map;
  }

  //   private static void printmap(HashMap<Integer, String> map) {
  ////    System.err.println("01234567890123456789012345678901234567890");
  ////    System.err.println("   012345678901234567890      12345678901234567890");
  ////    System.err.println("GACTTTGAGAGGGANAAAGTTATG      AACATTTATN");
  //    for (Map.Entry<Integer, String> e : map.entrySet()) {
  //     System.err.println(getPos(e.getKey()) + e.getValue());
  //     System.err.println("key = " + e.getKey());
  //    }
  //  }

  //  private static String getPos(int key) {
  //    final StringBuilder sb = new StringBuilder();
  //    for (int i = 0; i < key; ++i) {
  //      sb.append(' ');
  //    }
  //    return sb.toString();
  //  }

  static String unrollCigar(final SAMRecord rec) {
    final String cigar = rec.getCigarString();
    final String read = rec.getReadString().toUpperCase(Locale.getDefault());
    int n = 0;
    int zeroBasedreadPos = 0;
    final StringBuilder sb = new StringBuilder();
    for (int i = 0; i < cigar.length(); ++i) {
      final char c = cigar.charAt(i);
      if (Character.isDigit(c)) {
        n = 10 * n + c - '0';
      } else {
        for (int j = 0; j < n; ++j) {
          switch (c) {
            case SamUtils.CIGAR_SAME_OR_MISMATCH:
            case SamUtils.CIGAR_MISMATCH:
            case SamUtils.CIGAR_SAME:
            case SuperCigar.UNKNOWN_TEMPLATE:
            case SuperCigar.UNKNOWN_READ:
              sb.append(read.charAt(zeroBasedreadPos++));
              break;

            case SamUtils.CIGAR_INSERTION_INTO_REF:
              sb.append(Character.toLowerCase(read.charAt(zeroBasedreadPos++)));
              break;
            case SamUtils.CIGAR_DELETION_FROM_REF:
              sb.append("-");
              break;
            case SamUtils.CIGAR_GAP_IN_READ:
              sb.append(" ");
              break;
            case SamUtils.CIGAR_SOFT_CLIP:
              ++zeroBasedreadPos;
              break;
            case SamUtils.CIGAR_HARD_CLIP:
              break;
            default:
              throw new IllegalStateException(INVALID_CIGAR + cigar);
          }
        }
        n = 0;
      }
    }
    return sb.toString();
  }

}
