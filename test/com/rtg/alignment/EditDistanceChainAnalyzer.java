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
package com.rtg.alignment;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import com.rtg.mode.DnaUtils;
import com.rtg.ngs.NgsParams;
import com.rtg.ngs.NgsParamsBuilder;
import com.rtg.util.IntegerOrPercentage;
import com.rtg.util.MaxShiftFactor;
import com.rtg.util.StringUtils;
import com.rtg.util.cli.CFlags;
import com.rtg.util.diagnostic.Diagnostic;

import junit.framework.TestCase;

/**
 * This tests the whole alignment chain.
 * It is similar to prioritised edit distance, but tests the correctness of
 * ALL the aligners it is given, and measures their speeds.
 *
 * It is more than just a unit test (see its main method),
 * but it does include a few quick unit tests for regression purposes.
 *
 */
public class EditDistanceChainAnalyzer extends TestCase implements UnidirectionalEditDistance {
  /** Extra context bases printed either side of the template. */
  private static final int CONTEXT = 10;

  private final UnidirectionalEditDistance[] mEds;

  /** result from each aligner */
  private final int[][] mResult;

  /** nanoseconds taken by each aligner */
  private final long[] mTimeTaken;

  /** if non-null, records details about the alignment scores of each read */
  private boolean mReadDetails;

  /** if non-null, records details about each error that is detected */
  private boolean mErrorDetails;

  /** if non-null, log statistics about every aligner. */
  private boolean mVerboseStats;

  /** the maximum amount we try to move the start position */
  private int mShift;

  /** the name of the current read being aligned */
  private String mReadName;

  static final int RESULT_NULL = 0;
  static final int RESULT_MAX_CORRECT = 1;
  static final int RESULT_MAX_WRONG = 2;
  static final int RESULT_YES_TOO_LOW = 3;
  static final int RESULT_YES_CORRECT = 4;
  static final int RESULT_YES_TOO_HIGH = 5;
  static final int RESULT_EVENTS = 6;

  /** records counts of RESULT_* events for each aligner */
  private final int[][] mStats;

  /** total nanoseconds taken by each aligner */
  private final long[] mTotalTimeTaken;

  private final int[] mMaxIntActions;

  private long mTotalReadLen;
  private long mTotalMaxScore;

  private final int mMismatchPenalty;
  private final int mUnknownsPenalty;
  private final int mGapOpenPenalty;
  private final int mGapExtendPenalty;
  private final MaxShiftFactor mMaxShiftFactor;

  /**
   * Uses the aligners in the given order, but in self-testing mode.
   *
   * @param mismatchPenalty penalty score for a mismatch during alignment
   * @param gapOpenPenalty penalty score for a gap open during alignment
   * @param gapExtendPenalty penalty score for a gap extension during alignment
   * @param unknownsPenalty penalty score for an unknown nucleotide during alignment
   * @param editDistances the edit distances to iterate over (in order)
   * @param maxShiftFactor how far reads can shift when aligning
   */
  public EditDistanceChainAnalyzer(int mismatchPenalty, int gapOpenPenalty, int gapExtendPenalty, int unknownsPenalty, MaxShiftFactor maxShiftFactor, UnidirectionalEditDistance... editDistances) {
    mReadDetails = false;
    mErrorDetails = false;
    mEds = editDistances;
    mResult = new int[editDistances.length][];
    mTimeTaken = new long[editDistances.length];
    mTotalTimeTaken = new long[editDistances.length];

    mStats = new int[editDistances.length][];
    for (int i = 0; i < mStats.length; i++) {
      mStats[i] = new int[RESULT_EVENTS];
    }
    mMaxIntActions = new int[12];
    mMaxIntActions[ActionsHelper.ALIGNMENT_SCORE_INDEX] = Integer.MAX_VALUE;

    mMismatchPenalty = mismatchPenalty;
    mGapOpenPenalty = gapOpenPenalty;
    mGapExtendPenalty = gapExtendPenalty;
    mUnknownsPenalty = unknownsPenalty;
    mMaxShiftFactor = maxShiftFactor;
  }

  /**
   * Uses the following chain: no indels, lower-bound, hop-step-long, seeded, gotoh.
   *
   * @param params {@link NgsParams} for current run
   * @param rlen    Expected read lengths
   * @param maxShiftFactor controls aligner band widths
   */
  public EditDistanceChainAnalyzer(NgsParams params, int rlen, MaxShiftFactor maxShiftFactor) {
    this(params.substitutionPenalty(), params.gapOpenPenalty(), params.gapExtendPenalty(), params.unknownsPenalty(), maxShiftFactor,
        new NoIndelsEditDistance(params),
        new SingleIndelEditDistance(params, rlen),
//        new SingleIndelSeededEditDistance(params, rlen), XXX this is a broken unused aligner
        new LowerBoundEditDistance(EditDistanceFactory.calculateLowerBoundValue(rlen, maxShiftFactor), 1, params.unknownsPenalty()),
        new HopStepEditDistanceLong(params),
        new SeededAligner(params, false),
        new GotohEditDistance(params)
        );
  }

  /**
   * Uses gap open penalty of 1, read length 1000, and maximum score 5%.
   * @param params {@link NgsParams} for current run
   */
  public EditDistanceChainAnalyzer(NgsParams params) {
    this(params, 1000, new MaxShiftFactor(0.5));
  }

  /**
   * Print the results of each aligner for each read.
   *
   * @param eachRead true to turn on the printing.
   */
  public void printEachRead(final boolean eachRead) {
    mReadDetails = eachRead;
  }

  /**
   * Print details about each error by an aligner.
   *
   * @param eachError true to turn on the error details.
   */
  public void printEachError(final boolean eachError) {
    mErrorDetails = eachError;
  }

  /**
   * Try all start positions for each alignment, plus/minus <code>shift</code>.
   * @param shift the maximum shift to try.
   */
  public void setShift(final int shift) {
    mShift = shift;
  }

  /**
   * Print statistics for ALL aligners into the developer log.
   * @param verbose true means lots of output
   */
  public void setVerboseStats(final boolean verbose) {
    mVerboseStats = verbose;
  }

  @Override
  public int[] calculateEditDistance(byte[] read, int rlen, byte[] template, int zeroBasedStart, int maxScore, int maxShift, boolean cgLeft) {
    if (zeroBasedStart < 0 || rlen < 0 || zeroBasedStart >= template.length) {
      System.err.println("PROBLEM: read is off template: rlen=" + rlen + " start=" + zeroBasedStart
          + " len=" + template.length + " maxScore=" + maxScore + " read="
          + DnaUtils.bytesToSequenceIncCG(read));
    }
    int[] result = null;
    for (int ed = 0; ed < mEds.length; ed++) {
      final long starttime = System.nanoTime();
      mResult[ed] = mEds[ed].calculateEditDistance(read, rlen, template, zeroBasedStart, maxScore, maxShift, cgLeft);
      mTimeTaken[ed] = System.nanoTime() - starttime;
      if (result == null && mResult[ed] != null) {
        result = mResult[ed];
      }
    }
    analyzeResults(read, rlen, template, zeroBasedStart, maxScore);
    return result == null ? mMaxIntActions : result;
  }

  @Override
  public int[] calculateEditDistanceFixedBoth(byte[] read, int readStartPos, int readEndPos, byte[] template,
      int templateStartPos, int templateEndPos, int maxScore, int maxShift) {
    int[] result = null;
    for (int ed = 0; ed < mEds.length; ed++) {
      final long starttime = System.nanoTime();
      mResult[ed] = mEds[ed].calculateEditDistanceFixedBoth(read, readStartPos, readEndPos, template, templateStartPos, templateEndPos, maxScore, maxShift);
      mTimeTaken[ed] = System.nanoTime() - starttime;
      if (result == null && mResult[ed] != null) {
        result = mResult[ed];
      }
    }
    analyzeResults(read, readEndPos - readStartPos, template, templateStartPos, maxScore);
    return result == null ? mMaxIntActions : result;
  }

  @Override
  public int[] calculateEditDistanceFixedEnd(byte[] read, int readStartPos, int readEndPos, byte[] template,
      int templateExpectedStartPos, int templateEndPos, int maxScore, int maxShift) {
    int[] result = null;
    for (int ed = 0; ed < mEds.length; ed++) {
      final long starttime = System.nanoTime();
      mResult[ed] = mEds[ed].calculateEditDistanceFixedEnd(read, readStartPos, readEndPos, template, templateExpectedStartPos, templateEndPos,
          maxScore, maxShift);
      mTimeTaken[ed] = System.nanoTime() - starttime;
      if (result == null && mResult[ed] != null) {
        result = mResult[ed];
      }
    }
    analyzeResults(read, readEndPos - readStartPos, template, templateExpectedStartPos, maxScore);
    return result == null ? mMaxIntActions : result;
  }

  @Override
  public int[] calculateEditDistanceFixedStart(byte[] read, int readStartPos, int readEndPos, byte[] template,
      int templateStartPos, int maxScore, int maxShift) {
    int[] result = null;
    for (int ed = 0; ed < mEds.length; ed++) {
      final long starttime = System.nanoTime();
      mResult[ed] = mEds[ed].calculateEditDistanceFixedStart(read, readStartPos, readEndPos, template, templateStartPos, maxScore, maxShift);
      mTimeTaken[ed] = System.nanoTime() - starttime;
      if (result == null && mResult[ed] != null) {
        result = mResult[ed];
      }
    }
    analyzeResults(read, readEndPos - readStartPos, template, templateStartPos, maxScore);
    return result == null ? mMaxIntActions : result;
  }

  // analyze and record the results of this read.
  private void analyzeResults(byte[] read, int rlen, byte[] template, int zeroBasedStart, int maxScore) {
    final ActionsValidator av = new ActionsValidator(mGapOpenPenalty, mGapExtendPenalty, mMismatchPenalty, mUnknownsPenalty);
    final StringBuilder scoresLine = new StringBuilder();
    final StringBuilder startsLine = new StringBuilder();
    final StringBuilder timesLine = new StringBuilder();
    scoresLine.append(mReadName).append("\tscores");
    startsLine.append(mReadName).append("\tstarts");
    timesLine.append(mReadName).append("\tusecs");
    mTotalReadLen += rlen;
    mTotalMaxScore += maxScore;
    // now compare and validate the results
    final int[] last = mResult[mResult.length - 1]; //the last aligner (usually gotoh)'s actions array result

    final int lastScore = last == null ? 1000000 : last[ActionsHelper.ALIGNMENT_SCORE_INDEX]; //last aligner's score for the alignment
    boolean error = false;
    for (int ed = 0; ed < mEds.length; ed++) {
      mTotalTimeTaken[ed] += mTimeTaken[ed];
      timesLine.append("\t").append(mTimeTaken[ed] / 1000);
      final int[] actions = mResult[ed];
      if (actions == null) {
        mStats[ed][RESULT_NULL]++;
        scoresLine.append("\t-");
        startsLine.append("\t-");
      } else {
        // validate the results
        try {
          if (!av.isValid(actions, read, rlen, template, maxScore)) {
            error = true;
            final String msg = av.getErrorDetails(actions, read, rlen, template, zeroBasedStart);
            System.err.println("EditDistanceChainTest validation problem: " + getName(ed));
            System.err.println(msg);
          }
        } catch (final RuntimeException e) {
          error = true;
          System.err.println("EditDistanceChainTest: Problem even validating this pair of reads");
          System.err.println(" read:  " + DnaUtils.bytesToSequenceIncCG(read, 0, rlen));
          System.err.println(" tmpl:  " + DnaUtils.bytesToSequenceIncCG(template, zeroBasedStart, rlen));
          //e.printStackTrace();
        }
        final int score = actions[ActionsHelper.ALIGNMENT_SCORE_INDEX];
        if (score == Integer.MAX_VALUE) {
          if (lastScore == score) {
            mStats[ed][RESULT_MAX_CORRECT]++;
          } else {
            mStats[ed][RESULT_MAX_WRONG]++;
            error = true;
          }
          scoresLine.append("\tmax");
          startsLine.append("\tmax");
        } else {
          if (score == lastScore) {
            mStats[ed][RESULT_YES_CORRECT]++;
          } else if (score < lastScore) {
            mStats[ed][RESULT_YES_TOO_LOW]++;
//            if (lastScore != Integer.MAX_VALUE) {
//              System.err.println("ed " + mEds[ed].getClass().getName() + " : Lower than gotoh? g: " + lastScore + " s: " + score + " ts: " + zeroBasedStart + " r: " + DnaUtils.bytesToSequenceIncCG(read) + " t: " + DnaUtils.bytesToSequenceIncCG(template, zeroBasedStart, read.length));
//            }
            error = true;
          } else {
            mStats[ed][RESULT_YES_TOO_HIGH]++;
            error = true;
          }
          scoresLine.append("\t").append(actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
          startsLine.append("\t").append(actions[ActionsHelper.TEMPLATE_START_INDEX]);
        }
      }
    }
    if (mReadDetails || mErrorDetails && error) {
      System.out.println();
      System.out.println(scoresLine.toString() + (last == null ? "" : "\t0"));
      System.out.println(startsLine.toString() + "\t" + zeroBasedStart);
      System.out.println(timesLine.toString());
    }
    if (mErrorDetails && error) {
      final String space = StringUtils.getSpaceString(CONTEXT);
      System.out.println("\tread:\t" + space + DnaUtils.bytesToSequenceIncCG(read, 0, rlen));
      System.out.println("\ttmpl:\t" + DnaUtils.bytesToSequenceIncCG(template, zeroBasedStart - CONTEXT, rlen + 2 * CONTEXT));
      for (int ed = 0; ed < mEds.length; ed++) {
        final int[] actions = mResult[ed];
        if (actions != null && actions[ActionsHelper.ALIGNMENT_SCORE_INDEX] != Integer.MAX_VALUE) {
          final String name = getName(ed).substring(0, 6);
          System.out.println("\t" + name + ":\t" + space + ActionsHelper.toString(actions));
        }
      }
    }
  }

  /**
   * @param ed the number of an aligner, from 0 up to <code>getNumAligners() - 1</code>
   * @return the readable name of that aligner
   */
  public String getName(final int ed) {
    return mEds[ed].getClass().getSimpleName();
  }

  /**
   * Print statistics to developer log.
   * If <code>allDetails</code> is false, just the summary table will be logged.
   */
  @Override
  public void logStats() {
    int totalReads = 0;
    for (int event = 0; event < RESULT_EVENTS; event++) {
      totalReads += mStats[0][event];
    }
    final StringBuilder sb = new StringBuilder();
    sb.append("EditDistanceChainTest statistics,");
    sb.append(" readlen: ").append(mTotalReadLen / totalReads);
    sb.append(" maxScore: ").append(mTotalMaxScore / totalReads);
    sb.append(StringUtils.LS);
    sb.append("ED#\tNull#\tMaxOk\tMaxBad\tTooLow\tCorrect\tTooHigh\tTotal\tAvUsecs").append(StringUtils.LS);
    for (int ed = 0; ed < mEds.length; ed++) {
      if (getName(ed).length() > 6) { //for obfuscater
        sb.append(getName(ed).substring(0, 6));
      } else {
        sb.append(getName(ed));
      }
      int total = 0;
      for (int event = 0; event < RESULT_EVENTS; event++) {
        sb.append("\t");
        sb.append(mStats[ed][event]);
        total += mStats[ed][event];
      }
      assert total == totalReads;
      sb.append("\t");
      sb.append(total);
      sb.append("\t");
      sb.append(mTotalTimeTaken[ed] / total / 1000);
      sb.append(StringUtils.LS);

//      final double timetaken = mTimeTaken[i] * 1024 / 1000000000.0;
//      final double speed = (int) (mCounts[i] / timetaken);
//      final int consumed = mCounts[i] - mNullReturned[i];
//      final double eatingspeed = (int) (consumed / timetaken);
//      sb.append("  [ " + (i + 1) + " ] calls= " + mCounts[i] + " nulls= " + mNullReturned[i] + " MaxInt= " + mMaxIntReturned[i] + " avg.score= "
//          + (mCounts[i] != 0 ? "" + (mTotalScore[i] / mCounts[i]) : "-") + " time(s)= " + timetaken + " speed(per s)=" + speed + " eating(per s)="
//          + eatingspeed + "   "
//          + mEds[i].getClass().getName());
    }
    Diagnostic.developerLog(sb.toString());
    if (mVerboseStats) {
      for (UnidirectionalEditDistance mEd : mEds) {
        mEd.logStats();
      }
    }
  }

  /**
   * Run this aligner chain on a text file of read-template pairs,
   * and report statistics about the correctness of all the aligners.
   *
   * Each line of the input file must have the format:
   * "# comment line..." or
   * "ED: read template zeroBasedStart"
   *
   * @param inputFile the file to read
   * @param maxScore maximum score allowed
   * @throws IOException if file is unreadable or in wrong format.
   */
  public void doFile(final File inputFile, final int maxScore) throws IOException {
    doFile(new BufferedReader(new FileReader(inputFile)), maxScore);
  }

  /**
   * Run this aligner chain on a text file of read-template pairs,
   * and report statistics about the correctness of all the aligners.
   *
   * @param input the input stream to read
   * @param maxScore maximum score allowed
   * @throws IOException if file is unreadable or in wrong format.
   */
  public void doFile(final BufferedReader input, final int maxScore) throws IOException {
    try {
      String line;
      int lineNum = 0;
      while ((line = input.readLine()) != null) {
        lineNum++;
        if (line.startsWith("#")) {
          continue;
        }
        final String[] words = StringUtils.split(line, ' ');
        if (words.length != 4 || !words[0].equals("ED:")) {
          throw new IOException("Input line format should be: 'ED: read template zeroBasedStart'");
        }
        final byte[] read = DnaUtils.encodeString(words[1]);
        final byte[] tmpl = DnaUtils.encodeString(words[2]);
        final int startPos = Integer.parseInt(words[3]);
        mReadName = String.valueOf(lineNum);
        for (int pos = startPos - mShift; pos <= startPos + mShift; pos++) {
          calculateEditDistance(read, read.length, tmpl, pos, maxScore, mMaxShiftFactor.calculateMaxShift(read.length), true);
        }
      }
    } finally {
      input.close();
    }
  }

  /**
   * get name test
   */
  public void testGetName() {
    assertEquals(LowerBoundEditDistance.class.getSimpleName(), getName(0));
  }

  private static final String ALL_STATS = "allStats";
  private static final String ALL_READS = "allReads";
  private static final String ERROR_DETAILS = "errorDetails";
  private static final String MAX_SHIFT = "maxShift";
  private static final String MAX_SCORE = "maxScore";

  /**
   * init flags
   * @param flags input flags
   */
  public static void initFlags(CFlags flags) {
    flags.registerRequired(Integer.class, "READLEN", "Expected length of each read");
    flags.registerRequired(File.class, "INPUT-FILE", "File containing read and template strings");
//    flags.registerRequired(File.class, "SAM-FILE", "A sam file");
//    flags.registerRequired(File.class, "TEMPLATE-SDF", "SDF of the template");
    flags.registerOptional(ALL_STATS, "show detailed statistics from all aligners");
    flags.registerOptional(ALL_READS, "show all aligner scores for each read");
    flags.registerOptional(ERROR_DETAILS, "print the read, template and alignments of all erroneous alignment results");
    flags.registerOptional(MAX_SHIFT, Integer.class, "INT", "try all start positions from start-delta...start+delta", 0);
    flags.registerOptional(MAX_SCORE, String.class, "INT_OR_PERCENT", "the maximum score that aligners are allowed to return", "5%");
//    flags.registerOptional("MAX_RECORDS", Integer.class, "INT", "Max records to test");
  }

  /**
   * main
   * @param args command line arguments
   * @throws IOException if error
   */
  public static void main(final String[] args) throws IOException {
    final CFlags flags = new CFlags();
    initFlags(flags);
    if (!flags.setFlags(args)) {
      System.exit(1);
    }
    final boolean allStats = flags.isSet(ALL_STATS);
    final boolean allReads = flags.isSet(ALL_READS);
    final boolean errorDetails = flags.isSet(ERROR_DETAILS);
    final int shift = (Integer) flags.getValue(MAX_SHIFT);
    final IntegerOrPercentage maxScore = new IntegerOrPercentage((String) flags.getValue(MAX_SCORE));
    final int readlen = (Integer) flags.getAnonymousValue(0);
    final File file = (File) flags.getAnonymousValue(1);
//    final File template = (File) flags.getAnonymousValue(2);
    final NgsParams params = new NgsParamsBuilder().gapOpenPenalty(19).gapExtendPenalty(1).substitutionPenalty(9).unknownsPenalty(9).create();
    final EditDistanceChainAnalyzer ed = new EditDistanceChainAnalyzer(params, readlen, new MaxShiftFactor(0.5));
    ed.printEachRead(allReads);
    ed.printEachError(errorDetails);
    ed.setShift(shift);
    ed.setVerboseStats(allStats);

//    int maxRecords = Integer.MAX_VALUE;
//    if (flags.isSet("MAX_RECORDS")) {
//      maxRecords = (Integer) flags.getValue("MAX_RECORDS");
//    }
//
//    ed.doSamFile(file, template, readlen, maxRecords);

    ed.doFile(file, maxScore.getValue(readlen));
    ed.logStats();
  }
}
