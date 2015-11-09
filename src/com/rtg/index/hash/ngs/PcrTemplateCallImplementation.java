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
import com.rtg.launcher.HashingRegion;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Computes minimum distance between candidate probes and template sequences..
 */
public class PcrTemplateCallImplementation extends IntegralAbstract implements TemplateCall {

  private NgsHashFunction mHashFunction;

  private final IndexSet mIndexes;

  private final byte[] mDistances;

  private Finder mHit;

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("PcrTemplateCallImplementation");
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mIndexes != null && mIndexes.size() > 0);
    Exam.assertTrue(mDistances != null && mDistances.length > 0);
    Exam.assertTrue(mHit != null);
    return true;
  }

  /**
   * @param indexes indexes selected by the hash function.
   * @param distances minimum match distance so far
   */
  public PcrTemplateCallImplementation(final IndexSet indexes, final byte[] distances) {
    mHit = makeHit();
    mDistances = distances;
    mIndexes = indexes;
    //System.err.println("numberReads=" + mNumberReads + " threshold=" + threshold);
  }


  private Finder makeHit() {
    return new Finder() {
      @Override
      public void found(final long readIdl) {
        final int readId = (int) readIdl;
        final int score = mHashFunction.fastScore(readId);
        //System.err.println(" readId=" + readId + " score=" + score);
        assert score > 0 && score <= Byte.MAX_VALUE : score;
        if (score < mDistances[readId]) {
          mDistances[readId] = (byte) score;
        }
      }
    };
  }

  @Override
  public PcrTemplateCallImplementation clone() throws CloneNotSupportedException {
    final PcrTemplateCallImplementation clone = (PcrTemplateCallImplementation) super.clone();
    clone.mHit = clone.makeHit();
    assert integrity();
    return clone;
  }

  @Override
  public TemplateCall threadClone(final HashingRegion region) {
    if (region != HashingRegion.NONE) {
      throw new UnsupportedOperationException();
    }
    try {
    return clone();
    } catch (final CloneNotSupportedException e) {
      throw new RuntimeException(e);
    }
  }

  @Override
  public void threadFinish() {
  }

  @Override
  public void set(final long name, final int length) {
    //do nothing
  }

  @Override
  public void setReverse(final boolean reverse) {
    // do nothing
  }

  @Override
  public boolean isReverse() {
    return false;
  }

  @Override
  public void setHashFunction(final NgsHashFunction hashFunction) {
    mHashFunction = hashFunction;
  }

  @Override
  public void done() {
    //do nothing
  }

  @Override
  public void endSequence() {
    //do nothing
  }

  @Override
  public void templateCall(final int endPosition, final long hash, final int index) throws IOException {
    //System.err.println("search index=" + index + " hash=" + Utils.toBits(hash));
    mIndexes.get(index).search(hash, mHit);
  }

  @Override
  public void logStatistics() {
    // do nothing
  }
}
