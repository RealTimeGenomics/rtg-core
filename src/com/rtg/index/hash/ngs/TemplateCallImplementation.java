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
package com.rtg.index.hash.ngs;


import java.io.IOException;

import com.rtg.index.Finder;
import com.rtg.index.IndexSet;
import com.rtg.index.IntSet;
import com.rtg.index.IntSetCaller;
import com.rtg.index.IntSetSingle;
import com.rtg.index.IntSetWindow;
import com.rtg.launcher.HashingRegion;
import com.rtg.mode.SequenceType;
import com.rtg.ngs.NgsParams;
import com.rtg.reader.CgUtils;
import com.rtg.reader.PrereadType;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Does the actions for each position when scanning template.
 */
public class TemplateCallImplementation extends IntegralAbstract implements TemplateCall {
  protected long mTemplateId;
  protected int mTemplateLength;
  protected boolean mRC;
  protected int mEndPosition;

  NgsHashFunction mHashFunction;

  private final long mNumberReads;

  private final int mErrorLimit;

  private final boolean mIsCGV1;
  private final boolean mIsProtein;

  private final IndexSet mIndexes;

  private OutputProcessor mOutputProcessor;

  private IntSet mIS = null;

  private Finder mHit;

  private long mHitStatistics = 0;

  private long mProcessStatistics = 0;

  private final int mIntSetWindow;

  /**
   * @param params other configuration parameters
   * @param maxId the highest <code>internalId</code> that should be expected
   * @param indexes indexes selected by the hash function.
   * @param out where to place results.
   */
  public TemplateCallImplementation(final NgsParams params, final long maxId, final IndexSet indexes, final OutputProcessor out) {
    //System.err.println("entering");
    mHit = makeHit();
    mNumberReads = maxId;
    mIndexes = indexes;
    mOutputProcessor = out;
    mErrorLimit = params.outputParams().errorLimit();
    mIsCGV1 = params.buildFirstParams().reader().getPrereadType() == PrereadType.CG && params.buildFirstParams().reader().minLength() == CgUtils.CG_RAW_READ_LENGTH;
    mIntSetWindow = params.intSetWindow();
    mIsProtein = params.searchParams() != null && params.searchParams().reader().type() == SequenceType.PROTEIN;
  }

  private Finder makeHit() {
    //System.err.println("make Finder:" + System.identityHashCode(finder) + " parent TemplateCallImplementation:" + System.identityHashCode(this));
    return new Finder() {
      @Override
      public boolean found(final long readId) throws IOException {
        //System.err.println("adding readId=" + readId);
        //System.err.println("readid=" + readId + " " + mIS.toString());
        ++mHitStatistics;
        mIS.add((int) readId);
        return true;
      }
    };
  }

  @Override
  public TemplateCall threadClone(final HashingRegion region) throws IOException {
    try {
      final TemplateCallImplementation clone = (TemplateCallImplementation) super.clone();
      clone.mHashFunction = null;
      clone.mHit = clone.makeHit();
      clone.mOutputProcessor = mOutputProcessor.threadClone(region);
      //assert clone.integrity();
      return clone;
    } catch (final CloneNotSupportedException e) {
      throw new RuntimeException(e);
    }
  }

  public void setOutputProcessor(final OutputProcessor out) {
    mOutputProcessor = out;
  }

  @Override
  public void threadFinish() throws IOException {
    mOutputProcessor.threadFinish();
  }

  @Override
  public TemplateCallImplementation clone() throws CloneNotSupportedException {
    final TemplateCallImplementation clone = (TemplateCallImplementation) super.clone();
    clone.mHashFunction = null;
    clone.mHit = clone.makeHit();
    assert clone.integrity();
    return clone;
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("TemplateCallImplementation");
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mTemplateId >= 0);
    Exam.assertTrue(mTemplateLength >= 0);
    Exam.assertTrue(mEndPosition >= 0);
    Exam.assertTrue(mNumberReads >= 0);
    Exam.assertTrue(mIndexes != null && mIndexes.size() > 0);
    Exam.assertTrue(mErrorLimit >= 0);
    Exam.assertTrue(mHit != null);
    return true;
  }

  @Override
  public void set(final long name, final int length) {
    mTemplateId = name;
    mTemplateLength = length;
    if (mIS == null) {
      mIS = makeIntSet();
    }
  }

  @Override
  public void setReverse(final boolean reverse) {
    mRC = reverse;
  }

  @Override
  public boolean isReverse() {
    return mRC;
  }

  @Override
  public void setHashFunction(final NgsHashFunction hashFunction) {
    mHashFunction = hashFunction;
    mIS = null;
  }

  private IntSet makeIntSet() {
    final IntSetCaller caller = readId -> {
      final int scoreIndel = mHashFunction.indelScore(readId);
      //System.err.println("call score=" + score + " indelScore=" + scoreIndel + " errorLimit=" + mErrorLimit);
      if (scoreIndel <= mErrorLimit) {
        final int tzero = mIsProtein ? mEndPosition - mHashFunction.readLength() + 1 : Math.max(0, mEndPosition - mHashFunction.readLength() + 1);
        final String frame;
        if (mIsCGV1) {
          final boolean second = (readId & 1) == 1;
          frame = second == mRC ? "F" : "R";
        } else {
          frame = mRC ? "R" : "F";
        }
        ++mProcessStatistics;
        //System.err.println("calling process");
        final int score = mHashFunction.fastScore(readId);
        mOutputProcessor.process(mTemplateId, frame, readId, tzero, score, scoreIndel);
      }
    };
    int maxHashCount = mIndexes.get(0).maxHashCount();
    for (int i = 1; i < mIndexes.size(); ++i) {
      maxHashCount = Math.max(maxHashCount, mIndexes.get(i).maxHashCount());
    }
    //System.err.println("windows=" + mHashFunction.numberWindows() + " maxHashCount=" + mmaxHashCount);
    assert maxHashCount > 0;
    final IntSet intSet;
    if (mIntSetWindow == 1) {
      intSet = new IntSetSingle(mNumberReads, mHashFunction.numberWindows() * maxHashCount, caller);
    } else {
      //System.err.println("setting window");
      intSet = new IntSetWindow(mNumberReads, mHashFunction.numberWindows() * maxHashCount, mIntSetWindow, caller);
    }
    //System.err.println("make IntSet:" + System.identityHashCode(intSet) + " parent TemplateCallImplementation:" + System.identityHashCode(this));
    return intSet;
  }

  @Override
  public void done() throws IOException {
    //System.err.println("iterating over:" + System.identityHashCode(mIS) + " parent:" + System.identityHashCode(this));
    //System.err.println("done rc=" + mRC + " templateLength=" + mTemplateLength + " endPosition=" + mEndPosition);
//    mIndexes[0].dumpValues(System.err);
//    mIndexes[1].dumpValues(System.err);
//    mIndexes[2].dumpValues(System.err);
//    mIndexes[3].dumpValues(System.err);
//    mIndexes[4].dumpValues(System.err);
//    mIndexes[5].dumpValues(System.err);
    mIS.iterateClear();
  }

  /**
   * Flush any remaining values at an end of sequence.
   * @throws IOException If an error occurs
   */
  @Override
  public void endSequence() throws IOException {
    mIS.iterateClearAll();
  }

  @Override
  public void templateCall(final int endPosition, final long hash, final int index) throws IOException {
    mEndPosition = endPosition;
    //System.err.println("search index=" + index + " hash=" + hash + " ");

    mIndexes.get(index).search(hash, mHit);
  }

  @Override
  public void logStatistics() {
    Diagnostic.userLog("TemplateCall Hits " + mHitStatistics);
    Diagnostic.userLog("TemplateCall Processed " + mProcessStatistics);
  }
}

