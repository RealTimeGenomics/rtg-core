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


import java.io.ByteArrayOutputStream;

import com.rtg.alignment.ActionsHelper;
import com.rtg.variant.MachineErrorChooserInterface;
import com.rtg.variant.VariantAlignmentRecord;
import com.rtg.variant.VariantParams;
import com.rtg.variant.match.AlignmentMatch;
import com.rtg.variant.util.VariantUtils;

/**
 * This class to convert to and from CIGAR format, not remotely thread safe
 */
public final class CigarFormatter {

  private CigarFormatter() { }

  private static final int NONE = -1;
  private static final int MATCH_MISMATCH_LEGACY = 0;
  private static final int INS = 1;
  private static final int DEL = 2;
  private static final int SKIP = 3;
  private static final int SOFTCLIP = 4;
  private static final int MATCH = 7;
  private static final int MISMATCH = 8;

  private static final char[] CIGAR_CODES = SamUtils.getCigarCodes();

  private static void updateCigar(StringBuilder sb, int state, int count, boolean reverseComplement) {
    if (count > 0) {
      if (reverseComplement) {
        sb.insert(0, CIGAR_CODES[state]);
        sb.insert(0, count);
      } else {
        sb.append(count);
        sb.append(CIGAR_CODES[state]);
      }
    }
  }

  /**
   * Converts an array of actions and a read array into a cigar
   * Also soft-clips cigars if they extend past start/end of template.
   *
   *
   * @param actions the alignment actions
   * @param rc is the read reverse complement
   * @param templateLength length of the template sequence
   * @param legacy true if you want to use legacy cigar encoding (M for mismatch and match instead of X or =)
   * @param leftArm true for left arm, false for right arm (used only for CG Gotoh alignments)
   * @return a cigar string
   */
  public static String actionsToCigar(int[] actions, boolean rc, int templateLength, boolean legacy, boolean leftArm) {
    if (ActionsHelper.isCg(actions)) {
      return actionsToCigarCg(actions, rc, templateLength, legacy, leftArm);
    }
    final StringBuilder cigar = new StringBuilder();
    final ActionsHelper.CommandIterator iter = ActionsHelper.iterator(actions);
    int current = NONE;
    int previous = NONE;
    int count = 0;
    int tempPos = rc ? ActionsHelper.zeroBasedTemplateEndPos(actions) - 1 : ActionsHelper.zeroBasedTemplateStart(actions);
    final int direction = rc ? -1 : 1;

    while (iter.hasNext()) {
      final int action = iter.next();
      switch (action) {
        case ActionsHelper.SAME:
          current = tempPos < 0 || tempPos >= templateLength ? SOFTCLIP : MATCH;
          tempPos += direction;
          break;
        case ActionsHelper.INSERTION_INTO_REFERENCE:
          current = INS;
          break;
        case ActionsHelper.DELETION_FROM_REFERENCE:
          current = tempPos < 0 || tempPos >= templateLength ? SOFTCLIP : DEL;
          tempPos += direction;
          break;
        case ActionsHelper.MISMATCH:
        case ActionsHelper.UNKNOWN_TEMPLATE:
        case ActionsHelper.UNKNOWN_READ:
          current = tempPos < 0 || tempPos >= templateLength ? SOFTCLIP : MISMATCH;
          tempPos += direction;
          break;
        case ActionsHelper.SOFT_CLIP:
          current = SOFTCLIP;
          break;
        case ActionsHelper.NOOP:
          count--; //'continue' the previous action, but this wasn't actually a previous action.
          break;
        default:
          throw new IllegalArgumentException();
      }
      if (legacy && (current == MISMATCH || current == MATCH)) {
        current = MATCH_MISMATCH_LEGACY;
      }
      if (previous == NONE) {
        previous = current;
      }
      if (current == previous) {
        count++;
      } else {
        updateCigar(cigar, previous, count, rc);
        previous = current;
        count = 1;
      }
    }
    updateCigar(cigar, current, count, rc);
    return cigar.toString();
  }

  private static String actionsToCigarCg(int[] actions, boolean rc, int templateLength, boolean legacy, boolean leftArm) {
    final StringBuilder cigar = new StringBuilder();
    // We always process actions from the overlap end first, so that the commands AFTER the overlap
    // (towards the middle of the read) are ignored.
    final ActionsHelper.CommandIterator iter = leftArm ? ActionsHelper.iterator(actions) : ActionsHelper.iteratorReverse(actions);
    final boolean revCigar = leftArm == rc;  // very tricky -- one needs to think about all 4 cases here.
    // NOTE: the algorithm here is to process ALL the actions, to construct an ideal cigar and to find out
    //   how wide a portion of the template is spanned by the whole read (ignoring overhang).
    //   Then we postprocess the cigar to add soft clipping on the ends if necessary.
    int current = NONE;
    int previous = NONE;
    int count = 0;
    int cgOverlapSize = 0;

    int tempWidth = leftArm && !rc || !leftArm && rc ? actions[ActionsHelper.TEMPLATE_START_INDEX] : ActionsHelper.zeroBasedTemplateEndPos(actions) - 1;
    final int direction = leftArm && !rc || !leftArm && rc ? 1 : -1;

    while (iter.hasNext()) {
      final int action = iter.next();
      if (cgOverlapSize > 0 && action != ActionsHelper.CG_OVERLAP_IN_READ) {
        // start/continue skipping over cgOverlapSize nucleotides on the TEMPLATE.
        if (action != ActionsHelper.INSERTION_INTO_REFERENCE) {
          cgOverlapSize--;
          tempWidth += direction;
        }
        continue; // skip this command
      }
      switch (action) {
        case ActionsHelper.SAME:
          current = tempWidth < 0 || tempWidth >= templateLength ? SOFTCLIP : MATCH;
          tempWidth += direction;
          break;
        case ActionsHelper.CG_OVERLAP_IN_READ:
          tempWidth -= direction;
          count--; // preemptive attack against the increment down below, so that the cg overlap does not add to the cigar.
          // SAM Cigars do not handle CG overlap regions, so instead we will skip the appropriate number of actions after the overlap region.
          // Eg. If the actions are: ===OOSS===, we will skip the SS!
          cgOverlapSize++;
          break;
        case ActionsHelper.INSERTION_INTO_REFERENCE:
          current = INS;
          break;
        case ActionsHelper.CG_GAP_IN_READ:
          current = tempWidth < 0 || tempWidth >= templateLength ? SOFTCLIP : SKIP;
          tempWidth += direction;
          break;
        case ActionsHelper.DELETION_FROM_REFERENCE:
          current = tempWidth < 0 || tempWidth >= templateLength ? SOFTCLIP : DEL;
          tempWidth += direction;
          break;
        case ActionsHelper.MISMATCH:
        case ActionsHelper.UNKNOWN_TEMPLATE:
        case ActionsHelper.UNKNOWN_READ:
          current = tempWidth < 0 || tempWidth >= templateLength ? SOFTCLIP : MISMATCH;
          tempWidth += direction;
          break;
        case ActionsHelper.SOFT_CLIP:
          current = SOFTCLIP;
          break;
        case ActionsHelper.NOOP:
          count--;  //'continue' the previous action, but this wasn't actually a previous action.
          break;
        default:
          throw new IllegalArgumentException();
      }
      if (legacy && (current == MISMATCH || current == MATCH)) {
        current = MATCH_MISMATCH_LEGACY;
      }
      if (previous == NONE) {
        previous = current;
      }
      if (current == previous) {
        count++;
      } else if (previous == INS && current == DEL || previous == DEL && current == INS) {
        // adjacent Insert and Delete, so change 1D1I (or 1I1D) to 1X.
        updateCigar(cigar, previous, count - 1, revCigar);
        updateCigar(cigar, MISMATCH, 1, revCigar);
        previous = current;
        count = 0;
      } else {
        updateCigar(cigar, previous, count, revCigar);
        previous = current;
        count = 1;
      }
    }
    updateCigar(cigar, current, count, revCigar);
    return cigar.toString();
  }

  /**
   * Extract a match from a read which corresponds to a specified subsequence of the reference.
   * Care is needed to deal correctly with Ns and with partial overlaps by the read to the specified
   * subsequence. This can lead to all possibilities for open ended matches. In some cases a null may be returned
   * (for example when the subsequence corresponds to a mixture of Ns and other nucleotides).
   * The handling of inserts and deletions is particularly tricky here. An <a href="doc-files/indels.jpg">illustrative diagram</a> is available.
   * @param alignmentRecord the original sam record for the read.
   * @param chooser Machine error chooser (if null, no quality correction is done)
   * @param start position of the read on the template (0 based).
   * @param end position of the read on the template (0 based, exclusive).
   * @param params parameters from command line.
   * @return a match corresponding to the selected subset of the read (may be null).
   */
  public static AlignmentMatch cigarSubsequence(final VariantAlignmentRecord alignmentRecord, final MachineErrorChooserInterface chooser, final int start, final int end, final VariantParams params) {
    final String cigar = alignmentRecord.getCigar();
    final byte[] read = alignmentRecord.getRead();
    final byte[] qual = alignmentRecord.getRecalibratedQuality();
    if (read.length == 0) {
      return null; // For records that are non-primary and have no read or quality data stored with them
    }
    //position in reference (0 based)
    int refPos = alignmentRecord.getStart();
    //position in read (0 based)
    int readPos = 0;
    int n = 0;
    int i = 0;
    //move along cigar until reach start of selected subregion
    while (true) {
      if (i >= cigar.length()) {
        //ran out of room nothing to report
        return null;
      }
      final char c = cigar.charAt(i);
      if (Character.isDigit(c)) {
        n = 10 * n + c - '0';
      } else {
        while ((n > 0) && (refPos < start)) {
          switch (c) {
            case SamUtils.CIGAR_SAME_OR_MISMATCH:
            case SamUtils.CIGAR_SAME:
            case SamUtils.CIGAR_MISMATCH:
            case 'P':
              refPos++;
              readPos++;
              break;
            case SamUtils.CIGAR_INSERTION_INTO_REF:
            case SamUtils.CIGAR_SOFT_CLIP: // soft-clipping bases in read ignored for position
              readPos++;
              break;
            case SamUtils.CIGAR_GAP_IN_READ:
            case SamUtils.CIGAR_DELETION_FROM_REF:
              refPos++;
              break;
            case SamUtils.CIGAR_HARD_CLIP:
              break;
            default:
              throw new IllegalArgumentException(cigar);
          }
          n--;
        } //end while
        if ((n > 0) && (refPos >= start)) { // have reached selected region
          break;
        }
        assert n == 0;
      }
      i++;
    }
    assert n != 0;
    final int startInRead = readPos;
    //scan through the subregion setting flags about what has been seen
    //and building up the string for the subregion
    final StringBuilder sb = new StringBuilder();
    //final StringBuilder qb = qual.length == 0 ? null : new StringBuilder();
    final ByteArrayOutputStream qb = qual.length == 0 ? null : new ByteArrayOutputStream();
    boolean leftN = refPos > start;
    boolean validNt = false;
    boolean rightN = false;
    while (true) {
      if (i >= cigar.length()) {
        rightN = refPos < end;
        break;
      }
      final char c = cigar.charAt(i);
      if (Character.isDigit(c)) {
        n = 10 * n + c - '0';
      } else {
        while ((n > 0) && (c == SamUtils.CIGAR_INSERTION_INTO_REF ? refPos <= end : (refPos < end))) {
          switch (c) {
            case SamUtils.CIGAR_SAME_OR_MISMATCH:
            case SamUtils.CIGAR_SAME:
            case SamUtils.CIGAR_MISMATCH:
            case 'P'://never seen RTG doesn't generate it
              if (rightN) {
                //an invalid mixed case
                return null;
              }
              refPos++;
              //this duplicates the 'I' case below
              sb.append((char) read[readPos]);
              if (qb != null) {
                qb.write(qual[readPos]);
              }
              readPos++;
              //valid
              validNt = true;
              break;
            case SamUtils.CIGAR_INSERTION_INTO_REF:
              if (rightN) {
                //an invalid mixed case
                return null;
              }
              sb.append((char) read[readPos]);
              if (qb != null) {
                qb.write(qual[readPos]);
              }
              readPos++;
              //valid
              validNt = true;
              break;
            case SamUtils.CIGAR_SOFT_CLIP: // soft-clipping bases in read ignored for position
              readPos++;
              break;
            case SamUtils.CIGAR_GAP_IN_READ:
              refPos++;
              if (validNt) {
                rightN = true;
              } else {
                leftN = true;
              }
              break;
            case SamUtils.CIGAR_DELETION_FROM_REF:
              refPos++;
              break;
            case SamUtils.CIGAR_HARD_CLIP:
              break;
            default:
              throw new IllegalArgumentException("Unsupported cigar operation: " + c); //(c + ":" + cigar);
          }
          n--;
        }
        if ((n > 0) && (refPos >= end)) { // have reached selected region
          break;
        }
        assert n == 0 : n;
      }
      i++;
    }
    final int endReadPos = readPos;
    //distinguish case when insert at start and start and end are the same from there being no insert at start
    if (readPos == 0 && refPos == start) {
      return null;
    }
    if (leftN && !validNt && !rightN) {
      //all Ns - no useful match can be returned.
      return null;
    }
    final byte[] quality = qb == null ? null : qb.toByteArray();
    final String rs = sb.toString();
    final AlignmentMatch match = new AlignmentMatch(alignmentRecord, chooser, rs, quality, params.qDefault(), 0, rs.length(), VariantUtils.readScoreFromAlignmentRecord(alignmentRecord, params), !leftN, !rightN);
    match.setBasesLeftOfMatch(startInRead);
    match.setBasesRightOfMatch(alignmentRecord.getRead().length - endReadPos);
    return match;
  }
}
