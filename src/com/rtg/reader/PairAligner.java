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
class PairAligner {
  private final EditDistance mEd;
  private final int mMinOverlap;
  private final int mMinIdentity;
  private final ReadTrimmer mLengthTrimmer;
  private final PairAlignmentStats mStats = new PairAlignmentStats();
  private final int mR1ProbeLength;
  private final int mR2ProbeLength;
  private final boolean mVerbose;
  private final boolean mTrimMid;

  PairAligner(EditDistance ed, int minOverlap, int minIdentity, int r1ProbeLength, int r2ProbeLength, Integer minLength, boolean trimMid, boolean verbose) {
    mEd = ed;
    mMinOverlap = minOverlap;
    mMinIdentity = minIdentity;
    mR1ProbeLength = r1ProbeLength;
    mR2ProbeLength = r2ProbeLength;
    mLengthTrimmer = minLength > 0 ? new MinLengthReadTrimmer(minLength) : new NullReadTrimmer();
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

  void alignReads(FastqSequence template, FastqSequence read) {
    // Start talking in alignment speak.
    // Template is the R1
    // Read is the R2
    mStats.mTotal++;
    final int maxLength = Math.max(template.length(), read.length());
    final int[] actions = mEd.calculateEditDistance(read.getBases(), read.length(), template.getBases(), 0, false, Integer.MAX_VALUE, maxLength + 1, false);
    if (actions == null) {
      mStats.mNoAlignment++;
    } else {

      int mcount = 0;
      int templatePos = actions[ActionsHelper.TEMPLATE_START_INDEX];
      int readPos = 0;
      int lastReadBaseWithinTemplate = 0;
      int firstReadBaseWithinTemplate = -1;
      int lastTemplateBaseWithinRead = 0;
      int firstTemplateBaseWithinRead = -1;
      final ActionsHelper.CommandIterator ci = new ActionsHelper.CommandIterator();
      ci.setActions(actions);
      while (ci.hasNext()) {
        switch (ci.next()) {
          case SAME:
            if (templatePos >= 0 && templatePos < template.length()
              && readPos >= 0 && readPos < read.length()) {
              lastReadBaseWithinTemplate = readPos;
              lastTemplateBaseWithinRead = templatePos;
              firstReadBaseWithinTemplate = firstReadBaseWithinTemplate == -1 ? readPos : firstReadBaseWithinTemplate;
              firstTemplateBaseWithinRead = firstTemplateBaseWithinRead == -1 ? templatePos : firstTemplateBaseWithinRead;
            }
            mcount++;
            templatePos++;
            readPos++;
            break;
          case MISMATCH:
          case UNKNOWN_TEMPLATE:
          case UNKNOWN_READ:
            templatePos++;
            readPos++;
            break;
          case INSERTION_INTO_REFERENCE:
            readPos++;
            break;
          case DELETION_FROM_REFERENCE:
            templatePos++;
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
        doTrimAndRecordStats(template, read, actions, firstTemplateBaseWithinRead, lastTemplateBaseWithinRead, firstReadBaseWithinTemplate, lastReadBaseWithinTemplate);
      } else {
        mStats.mPoorAlignment++;
      }
    }
  }

  private void doTrimAndRecordStats(FastqSequence r1, FastqSequence r2, int[] actions, int overlapStartWithinR1, int overlapEndWithinR1, int overlapStartWithinR2, int overlapEndWithinR2) {
    if (mVerbose) {
      dumpAlignment(" IN-", r1, r2, actions);
    }
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

  private int countLeadingInsertions(int[] actions) {
    final ActionsHelper.CommandIterator ci = new ActionsHelper.CommandIterator();
    ci.setActions(actions);
    int count = 0;
    while (ci.hasNext()) {
      switch (ci.next()) {
        case SAME:
        case MISMATCH:
        case UNKNOWN_TEMPLATE:
        case UNKNOWN_READ:
        case DELETION_FROM_REFERENCE:
          return count;
        case INSERTION_INTO_REFERENCE:
          count++;
          break;
        default:
          throw new RuntimeException("Unhandled action encountered in alignment: " + ActionsHelper.toString(actions));
      }
    }
    return count;
  }

  protected void dumpAlignment(String prefix, FastqSequence r1, FastqSequence r2, int[] actions) {
    dumpAlignment(prefix, r1, r2, actions, 0, 0);
  }

  protected void dumpAlignment(String prefix, FastqSequence r1, FastqSequence r2, int[] actions, int r1Pad, int r2Pad) {
    final int pos = actions[ActionsHelper.TEMPLATE_START_INDEX];
    final int pad = Math.max(0, -pos);
    System.err.println(prefix + "R1: " + StringUtils.spaces(r1Pad + pad + countLeadingInsertions(actions)) + DnaUtils.bytesToSequenceIncCG(r1.getBases()));
    System.err.println(prefix + "AL: " + StringUtils.spaces(pad + pos) + ActionsHelper.toString(actions));
    System.err.println(prefix + "R2: " + StringUtils.spaces(r2Pad + pad + pos) + DnaUtils.bytesToSequenceIncCG(r2.getBases()));
  }

  public PairAlignmentStats getStats() {
    return mStats;
  }
}
