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
  private final PairAlignmentStats mStats = new PairAlignmentStats();
  private final int mProbeLength;

  PairAligner(EditDistance ed, int minOverlap, int minIdentity, int probeLength) {
    mEd = ed;
    mMinOverlap = minOverlap;
    mMinIdentity = minIdentity;
    mProbeLength = probeLength;
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
//      System.err.println("firstTemplateBaseWithinRead=" + firstTemplateBaseWithinRead);
//      System.err.println("lastTemplateBaseWithinRead=" + lastTemplateBaseWithinRead);
//      System.err.println("firstReadBaseWithinTemplate=" + firstReadBaseWithinTemplate);
//      System.err.println("lastReadBaseWithinTemplate=" + lastReadBaseWithinTemplate);

      //dumpAlignment(template, read, actions);
      if (overlap >= Math.min(0.8 * maxLength, mMinOverlap) && ident >= mMinIdentity) {

//        dumpAlignment(template, read, actions);

        doTrimAndRecordStats(template, read, firstReadBaseWithinTemplate, lastTemplateBaseWithinRead, firstTemplateBaseWithinRead);
      } else {
        mStats.mPoorAlignment++;
      }
    }
  }

  private void doTrimAndRecordStats(FastqSequence r1, FastqSequence r2, int overlapStartWithinR2, int overlapEndWithinR1, int overlapStartWithinR1) {
    if (overlapStartWithinR2 > 0) {
      // Read-through has occurred, remove the r2 bases
      final int trimAmount = overlapStartWithinR2 + mProbeLength;
      r2.trim(new FirstBasesReadTrimmer(trimAmount));
      mStats.mReadThroughOnR2++;
    } else if (overlapStartWithinR1 < mProbeLength) {
      // Large overlap but R2 probably contains probe bases.
      final int trimAmount = mProbeLength - overlapStartWithinR1;
      r2.trim(new FirstBasesReadTrimmer(trimAmount));
      mStats.mMajorOverlap++;
    } else {
      // Regular read overlap -- perhaps remove bases off the longest side
      mStats.mPartialOverlap++;
    }
    if (overlapEndWithinR1 < r1.length() - 1) {
      // Trim readthrough on the r1.
      mStats.mReadThroughOnR1++;
      final int numBases = r1.length() - overlapEndWithinR1 - 1;
      r1.trim(new LastBasesReadTrimmer(numBases));
    }
  }

  protected void dumpAlignment(FastqSequence r1, FastqSequence r2, int[] actions) {
    final int pos = actions[ActionsHelper.TEMPLATE_START_INDEX];
    final int pad = Math.max(0, -pos);
    System.err.println(StringUtils.spaces(pad + pos) + ActionsHelper.toString(actions));
    System.err.println(StringUtils.spaces(pad) + DnaUtils.bytesToSequenceIncCG(r1.getBases()));
    System.err.println(StringUtils.spaces(pad + pos) + DnaUtils.bytesToSequenceIncCG(r2.getBases()));
  }

  public PairAlignmentStats getStats() {
    return mStats;
  }
}
