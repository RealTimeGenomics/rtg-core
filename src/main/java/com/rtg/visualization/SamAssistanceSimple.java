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

package com.rtg.visualization;

import com.rtg.launcher.globals.CoreGlobalFlags;
import com.rtg.launcher.globals.GlobalFlags;
import com.rtg.sam.SamUtils;
import com.rtg.util.StringUtils;

import htsjdk.samtools.SAMRecord;

/**
 */
public class SamAssistanceSimple implements SamAssistance {

  private static final String INVALID_CIGAR = "Invalid cigar: ";
  private static final boolean PRETTY_DELETES = GlobalFlags.getBooleanValue(CoreGlobalFlags.CP437_DELETES);

  @Override
  public String[] samToReads(final SAMRecord sam, final String template, byte[] templateBytes, final int readStart, final boolean displayDots, boolean displaySoftClip) {
    final String read = sam.getReadString();
    if ("*".equals(read)) {
      throw new IllegalStateException("No read sequence");
    }
    final String cigar = sam.getCigarString();
    return new String[] {cigarToReads(cigar, read, template, readStart, displayDots, displaySoftClip)};
  }

  /**
   * @param cigar string.
   * @param read nucleotides of read.
   * @param template nucleotides of template
   * @param readStart where the read starts on the template (in screen co-ordinates).
   * @param displayDots if lines match display dot
   * @param displaySoftClip if soft clipped bases should be displayed
   * @return return a string with a display of the current read.
   */
  String cigarToReads(String cigar, String read, String template, int readStart, boolean displayDots, boolean displaySoftClip) {
    //System.err.println("cigar=" + cigar + " read=" + read + " template=" + template + " readStart=" + readStart);
    final StringBuilder sb = new StringBuilder();
    int n = 0;
    int rPos = 0;
    //beware tPos is in screen co-ordinates
    int tPos = readStart;
    //assert template.charAt(tPos) != '_';
    boolean isFirstAction = true; //this is so we can deal properly with inserts at the start of reads... although it's sure to be broken for reads that start with inserts that aren't at the start of the template.
    for (int i = 0; i < cigar.length(); ++i) {
      final char c = cigar.charAt(i);
      if (Character.isDigit(c)) {
        n = 10 * n + c - '0';
      } else {
        assert n > 0 : cigar;
        //System.err.println("i=" + i + " n=" + n + " c=" + c);
        for (int j = 0; j < n; ) {
          if (tPos >= template.length()) {
            break;  //oh shiiiiiii
          }
          //System.err.println("  j=" + j + " rPos=" + rPos + " tPos=" + tPos + " te=" + template.charAt(tPos) + " sb=" + sb.toString());
          if (c == SamUtils.CIGAR_INSERTION_INTO_REF) {
            if (!isFirstAction) {
              tPos += n;
            } else {
              sb.append(StringUtils.spaces(readStart - n));
              isFirstAction = false;
            }
            // assert template.charAt(tPos) == '_' : tPos;
            for (int k = 0; k < n; ++k) {
              sb.append(read.charAt(rPos + k));
            }
            rPos += n;

            break;  //process ALL the inserts at the same time, so break out of this entire loop...
          }
          if (template.charAt(tPos) == DisplayHelper.INSERT_CHAR) {
            sb.append(DisplayHelper.INSERT_CHAR);
            ++tPos;
            continue;
          }
          ++j;

          switch (c) {
            case SamUtils.CIGAR_SAME:
            case SamUtils.CIGAR_SAME_OR_MISMATCH:
            case SamUtils.CIGAR_MISMATCH:
              if (isFirstAction) {
                sb.append(StringUtils.spaces(readStart));
                isFirstAction = false;
              }
              final char readChar = Character.toLowerCase(read.charAt(rPos));
              final char templateChar = Character.toLowerCase(template.charAt(tPos));
              if (readChar == templateChar && displayDots) {
                if (readChar == 'n') {
                  sb.append('N');
                } else {
                  sb.append('.');
                }
              } else {
                sb.append(read.charAt(rPos));
              }
              ++rPos;
              ++tPos;
              break;
            case SamUtils.CIGAR_DELETION_FROM_REF:
              if (isFirstAction) {
                sb.append(StringUtils.spaces(readStart));
                isFirstAction = false;
              }
              if (PRETTY_DELETES) {
                if (sb.charAt(sb.length() - 1) == '\u2551') {
                  sb.deleteCharAt(sb.length() - 1);
                  sb.append("\u251c\u2524");
                } else if (sb.charAt(sb.length() - 1) == '\u2524') {
                  sb.deleteCharAt(sb.length() - 1);
                  sb.append("\u254c\u2524");
                } else {
                  sb.append("\u2551");
                }
              } else {
                sb.append("-");
              }
              ++tPos;
              break;
            case SamUtils.CIGAR_GAP_IN_READ:
              if (isFirstAction) {
                sb.append(StringUtils.spaces(readStart));
                isFirstAction = false;
              }
              sb.append(" ");
              ++tPos;
              break;
            case SamUtils.CIGAR_SOFT_CLIP:
              if (displaySoftClip) {
                if (readStart - n + rPos >= 0) {
                  if (isFirstAction) {
                    sb.append(StringUtils.spaces(readStart - n + rPos));
                    isFirstAction = false;
                  }
                  sb.append(Character.toLowerCase(read.charAt(rPos)));
                }
              }
              ++rPos;
              //soft clip doesn't really count as a first action!
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
