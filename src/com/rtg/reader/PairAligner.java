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

package com.rtg.reader;

import static com.rtg.alignment.ActionsHelper.DELETION_FROM_REFERENCE;
import static com.rtg.alignment.ActionsHelper.INSERTION_INTO_REFERENCE;
import static com.rtg.alignment.ActionsHelper.MISMATCH;
import static com.rtg.alignment.ActionsHelper.SAME;
import static com.rtg.alignment.ActionsHelper.UNKNOWN_READ;
import static com.rtg.alignment.ActionsHelper.UNKNOWN_TEMPLATE;

import java.util.List;

import com.rtg.alignment.ActionsHelper;
import com.rtg.alignment.EditDistance;
import com.rtg.mode.DnaUtils;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;

/**
 */
public class PairAligner {

  /**
   * Select the behaviour when encountering a mismatch/indel within the overlap region
   */
  public enum MismatchType {
    /** Do nothing */
    NONE {
      @Override
      void apply(FastqSequence r1, FastqSequence r2, int[] actions) { }
    },

    /** Set the base quality of non-matching bases to zero */
    ZERO_PHRED {
      @Override
      void doMismatch(FastqSequence r1, FastqSequence r2, int r1Pos, int r2Pos) {
        r1.getQualities()[r1Pos] = 0;
        r2.getQualities()[r2Pos] = 0;
      }
      @Override
      void doIndel(FastqSequence r, int rPos) {
        r.getQualities()[rPos] = 0;
      }
    },

    /** Select the base with the highest reported quality, takes no action with indels */
    PICK_BEST {
      @Override
      void doMismatch(FastqSequence r1, FastqSequence r2, int r1Pos, int r2Pos) {
        final byte tq = r1.getQualities()[r1Pos];
        final byte rq = r2.getQualities()[r2Pos];
        if (tq > rq) {
          r2.getBases()[r2Pos] = r1.getBases()[r1Pos];
          r2.getQualities()[r2Pos] = tq;
        } else if (tq < rq) {
          r1.getBases()[r1Pos] = r2.getBases()[r2Pos];
          r1.getQualities()[r1Pos] = rq;
        }
      }
    };

    void apply(FastqSequence r1, FastqSequence r2, int[] actions) {
      int r1Pos = actions[ActionsHelper.TEMPLATE_START_INDEX];
      int r2Pos = 0;
      final ActionsHelper.CommandIterator ci = new ActionsHelper.CommandIterator();
      ci.setActions(actions);
      while (ci.hasNext()) {
        switch (ci.next()) {
          case SAME:
          case UNKNOWN_TEMPLATE:
          case UNKNOWN_READ:
            r1Pos++;
            r2Pos++;
            break;
          case MISMATCH:
            if (r1Pos >= 0 && r2Pos >= 0 && r1Pos < r1.length() && r2Pos < r2.length()) {
              doMismatch(r1, r2, r1Pos, r2Pos);
            }
            r1Pos++;
            r2Pos++;
            break;
          case INSERTION_INTO_REFERENCE:
            if (r1Pos >= 0 && r2Pos >= 0 && r1Pos < r1.length() && r2Pos < r2.length()) {
              doIndel(r2, r2Pos);
            }
            r2Pos++;
            break;
          case DELETION_FROM_REFERENCE:
            if (r1Pos >= 0 && r2Pos >= 0 && r1Pos < r1.length() && r2Pos < r2.length()) {
              doIndel(r1, r1Pos);
            }
            r1Pos++;
            break;
          default:
            throw new RuntimeException("Unhandled action encountered in alignment: " + ActionsHelper.toString(actions));
        }
      }
    }
    void doMismatch(FastqSequence r1, FastqSequence r2, int r1Pos, int r2Pos) { }
    void doIndel(FastqSequence r, int rPos) { }
  }

  private final EditDistance mEd;
  private final int mMinOverlap;
  private final int mMinIdentity;
  private final ReadTrimmer mLengthTrimmer;
  private final PairAlignmentStats mStats = new PairAlignmentStats();
  private final int mR1ProbeLength;
  private final int mR2ProbeLength;
  private final boolean mVerbose;
  private final boolean mTrimMid;
  private final MismatchType mMismatch;

  PairAligner(EditDistance ed, int minOverlap, int minIdentity, int r1ProbeLength, int r2ProbeLength, int minLength, boolean trimMid, MismatchType mismatch, boolean verbose) {
    mEd = ed;
    mMinOverlap = minOverlap;
    mMinIdentity = minIdentity;
    mR1ProbeLength = r1ProbeLength;
    mR2ProbeLength = r2ProbeLength;
    mLengthTrimmer = minLength > 0 ? new MinLengthReadTrimmer(minLength) : new NullReadTrimmer();
    mMismatch = mismatch;
    mTrimMid = trimMid;
    mVerbose = verbose;
  }

  void processReads(List<FastqPair> pairs) {
    Diagnostic.developerLog("Starting batch " + this);
    for (FastqPair pair : pairs) {
      processReads(pair);
    }
    Diagnostic.developerLog("Finished batch " + this);
  }

  void processReads(FastqPair pair) {
    // Reverse complement R2
    pair.r2().rc();

    // Look for high-identity alignment
    alignReads(pair.r1(), pair.r2());

    pair.r1().trim(mLengthTrimmer);
    pair.r2().trim(mLengthTrimmer);

    // Re-reverse R2 back to original order
    pair.r2().rc();
  }

  void alignReads(FastqSequence r1, FastqSequence r2) {
    // Start talking in alignment speak.
    // Template is the R1
    // Read is the R2
    mStats.mTotal++;
    final int maxLength = Math.max(r1.length(), r2.length());
    final int[] actions = mEd.calculateEditDistance(r2.getBases(), r2.length(), r1.getBases(), 0, false, Integer.MAX_VALUE, maxLength + 1, false);
    if (actions == null) {
      mStats.mNoAlignment++;
    } else {

      int mcount = 0;
      int r1Pos = actions[ActionsHelper.TEMPLATE_START_INDEX];
      int r2Pos = 0;
      int lastReadBaseWithinTemplate = 0;
      int firstReadBaseWithinTemplate = -1;
      int lastTemplateBaseWithinRead = 0;
      int firstTemplateBaseWithinRead = -1;
      final ActionsHelper.CommandIterator ci = new ActionsHelper.CommandIterator();
      ci.setActions(actions);
      while (ci.hasNext()) {
        switch (ci.next()) {
          case SAME:
            if (r1Pos >= 0 && r1Pos < r1.length()
              && r2Pos >= 0 && r2Pos < r2.length()) {
              lastReadBaseWithinTemplate = r2Pos;
              lastTemplateBaseWithinRead = r1Pos;
              firstReadBaseWithinTemplate = firstReadBaseWithinTemplate == -1 ? r2Pos : firstReadBaseWithinTemplate;
              firstTemplateBaseWithinRead = firstTemplateBaseWithinRead == -1 ? r1Pos : firstTemplateBaseWithinRead;
            }
            mcount++;
            r1Pos++;
            r2Pos++;
            break;
          case MISMATCH:
          case UNKNOWN_TEMPLATE:
          case UNKNOWN_READ:
            r1Pos++;
            r2Pos++;
            break;
          case INSERTION_INTO_REFERENCE:
            r2Pos++;
            break;
          case DELETION_FROM_REFERENCE:
            r1Pos++;
            break;
          default:
            throw new RuntimeException("Unhandled action encountered in alignment: " + ActionsHelper.toString(actions));
        }
      }

      final int overlap = Math.max(lastReadBaseWithinTemplate - firstReadBaseWithinTemplate,
        lastTemplateBaseWithinRead - firstTemplateBaseWithinRead) + 1; // End is not exclusive
      final int ident = overlap == 0 ? 0 : mcount * 100 / overlap;
      //System.out.println("Alignment score: " + score + " at " + pos + " with overlap " + overlap + " identity " + ident);
      //dumpAlignment(template, read, actions);
      if (overlap >= Math.min(0.8 * maxLength, mMinOverlap) && ident >= mMinIdentity) {
        doTrimAndRecordStats(r1, r2, actions, firstTemplateBaseWithinRead, lastTemplateBaseWithinRead, firstReadBaseWithinTemplate, lastReadBaseWithinTemplate);
      } else {
        mStats.mPoorAlignment++;
      }
    }
  }

  private void doTrimAndRecordStats(FastqSequence r1, FastqSequence r2, int[] actions, int overlapStartWithinR1, int overlapEndWithinR1, int overlapStartWithinR2, int overlapEndWithinR2) {
    if (mVerbose) {
      dumpAlignment(" IN-", r1, r2, actions);
    }

    mMismatch.apply(r1, r2, actions);

    final int r2Length = r2.length();
    mStats.mOverlapping++;
    mStats.mFragLengths.increment(overlapEndWithinR1 + (r2Length - overlapEndWithinR2));
    mStats.mOverlapDist.increment(overlapEndWithinR1 - overlapStartWithinR1);
    final int r2ProbeTrimAmount = mR1ProbeLength - overlapStartWithinR1;
    final int r1ProbeTrimAmount = overlapEndWithinR2 - (r2Length - mR2ProbeLength);
    if (overlapStartWithinR2 > 0) {
      // Read-through past start of R1 has occurred, remove the r2 bases
      final int trimAmount = overlapStartWithinR2 + mR1ProbeLength;
      r2.trim(new FirstBasesReadTrimmer(trimAmount));
      mStats.mR2ReadThrough++;
    } else if (r2ProbeTrimAmount > 0) {
      // Large overlap into R1 read start area, remove R2 bases.
      r2.trim(new FirstBasesReadTrimmer(r2ProbeTrimAmount));
      mStats.mR2ReadIntoR1Probe++;
    }
    if (overlapEndWithinR1 < r1.length() - 1) {
      // R1 read through past R2 start, trim the extra r1 bases.
      mStats.mR1ReadThrough++;
      final int numBases = r1.length() - overlapEndWithinR1 - 1 + mR2ProbeLength;
      r1.trim(new LastBasesReadTrimmer(numBases));
    } else if (r1ProbeTrimAmount > 0) {
      // Large overlap into R2 read start area, remove R1 bases.
      r1.trim(new LastBasesReadTrimmer(r1ProbeTrimAmount));
      mStats.mR1ReadIntoR2Probe++;
    }

    if (mTrimMid) {
      if (mVerbose) {
        dumpAlignment("INT-", r1, r2, actions, 0, r2Length - r2.length());
      }
      // Trim both sides to the midpoint of remaining overlap
      final int overlapR1Len = overlapEndWithinR1 - overlapStartWithinR1 + 1;
      final int newOverlapR2Len = overlapEndWithinR2 - Math.max(overlapStartWithinR2, r2Length - r2.length()) + 1;
      final int newOverlapR1Len = Math.min(overlapEndWithinR1, r1.length() - 1) - overlapStartWithinR1 + 1;
      //System.err.println("Overlap initial: R1=" + overlapR1Len + " R2=" + overlapR2Len + " remaining: R1=" + newOverlapR1Len + " R2=" + newOverlapR2Len);
      final int totalOverlap = (newOverlapR1Len + newOverlapR2Len) - overlapR1Len;
      int midTrim = totalOverlap / 2;
      if (midTrim > 0) {
        r1.trim(new LastBasesReadTrimmer(midTrim));
        if (midTrim * 2 < totalOverlap) {
          midTrim++;
        }
        r2.trim(new FirstBasesReadTrimmer(midTrim));
      }
    }
    if (mVerbose) {
      dumpAlignment("OUT-", r1, r2, actions, 0, r2Length - r2.length());
      System.err.println();
    }
  }

  /**
   * Convert a sequence and actions into a form suitable for display.
   * @param actions alignment actions
   * @param bases bases of sequence
   * @param ref true if this is the reference sequence
   * @param skipStart offset
   * @return read string
   */
  public static String bytesToSequence(int[] actions, byte[] bases, boolean ref, int skipStart) {
    final StringBuilder sb = new StringBuilder();
    int refPos = actions[ActionsHelper.TEMPLATE_START_INDEX];
    int readPos = 0;
    if (ref && refPos > 0) {
      sb.append(DnaUtils.bytesToSequenceIncCG(bases, 0, refPos));
    }

    readPos += ref ? 0 : -skipStart;
    refPos += ref ? -skipStart : 0;

    final ActionsHelper.CommandIterator ci = new ActionsHelper.CommandIterator();
    ci.setActions(actions);
    while (ci.hasNext()) {
      switch (ci.next()) {
        case SAME:
        case MISMATCH:
        case UNKNOWN_TEMPLATE:
        case UNKNOWN_READ:
          if (ref) {
            sb.append(base(bases, refPos));
          } else {
            sb.append(base(bases, readPos));
          }
          refPos++;
          readPos++;
          break;
        case DELETION_FROM_REFERENCE:
          sb.append(ref ? base(bases, refPos) : ' ');
          refPos++;
          break;
        case INSERTION_INTO_REFERENCE:
          sb.append(!ref ? base(bases, readPos) : ' ');
          readPos++;
          break;
        default:
          throw new RuntimeException("Unhandled action encountered in alignment: " + ActionsHelper.toString(actions));
      }
    }

    readPos += ref ? 0 : skipStart;
    refPos += ref ? skipStart : 0;

    if (ref && refPos < bases.length) {
      sb.append(DnaUtils.bytesToSequenceIncCG(bases, refPos, bases.length - refPos));
    } else if (!ref && readPos < bases.length) {
      throw new IllegalStateException("Not all read bases consumed by actions");
    }
    return sb.toString();
  }

  private static char base(final byte[] a, final int p) {
    return p >= 0 && p < a.length ? DnaUtils.getBase(a[p]) : ' ';
  }

  protected void dumpAlignment(String prefix, FastqSequence r1, FastqSequence r2, int[] actions) {
    dumpAlignment(prefix, r1, r2, actions, 0, 0);
  }

  protected void dumpAlignment(String prefix, FastqSequence r1, FastqSequence r2, int[] actions, int r1Pad, int r2Pad) {
    final int pos = actions[ActionsHelper.TEMPLATE_START_INDEX];
    final int extraRef = pos > 0 ? pos : 0;
    final String r1Str = bytesToSequence(actions, r1.getBases(), true, r1Pad);
    final String r2Str = bytesToSequence(actions, r2.getBases(), false, r2Pad);
    final String actionStr = ActionsHelper.toString(actions);
    System.err.println(prefix + "R1: " + r1Str);
    System.err.println(prefix + "AL: " + StringUtils.spaces(extraRef) + actionStr);
    System.err.println(prefix + "R2: " + StringUtils.spaces(extraRef) + r2Str);
  }

  public PairAlignmentStats getStats() {
    return mStats;
  }
}
