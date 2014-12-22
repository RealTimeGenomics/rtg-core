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
package com.rtg.protein;

import com.rtg.index.Index;
import com.rtg.index.IndexSet;
import com.rtg.index.hash.ngs.NgsHashFunction;
import com.rtg.index.hash.ngs.ProteinNgsHashLoop;
import com.rtg.index.hash.ngs.TemplateCall;
import com.rtg.launcher.HashingRegion;

import junit.framework.TestCase;

/**
 * Test class
 */
public class ReadLengthHashingStateTest extends TestCase {

  public void test() {
    final ProteinNgsHashLoop firstLoop = new ProteinNgsHashLoop(32, 32);
    final NgsHashFunction firstHash = new FakeHashFunction();
    final IndexSet firstIndex = new IndexSet(new Index[0]);
    final TemplateCall firstTci = new FakeTemplateCall();
    final ReadLengthHashingState rlhs = new ReadLengthHashingState(firstIndex, firstLoop, firstHash, firstTci);
    assertTrue(firstLoop == rlhs.getHashLoop());
    assertTrue(firstHash == rlhs.getHashFunction());
    assertTrue(firstIndex == rlhs.getIndexes());
    assertTrue(firstTci == rlhs.getTemplateCallImplementation());

  }
  private static class FakeTemplateCall implements TemplateCall {
    @Override
    public TemplateCall clone() throws CloneNotSupportedException {
      super.clone();
      throw new UnsupportedOperationException("Not supported yet.");
    }
    @Override
    public void done() {
      throw new UnsupportedOperationException("Not supported yet.");
    }
    @Override
    public void endSequence() {
      throw new UnsupportedOperationException("Not supported yet.");
    }
    @Override
    public boolean isReverse() {
      throw new UnsupportedOperationException("Not supported yet.");
    }
    @Override
    public void logStatistics() {
      throw new UnsupportedOperationException("Not supported yet.");
    }
    @Override
    public void set(long name, int length) {
      throw new UnsupportedOperationException("Not supported yet.");
    }
    @Override
    public void setHashFunction(NgsHashFunction hashFunction) {
      throw new UnsupportedOperationException("Not supported yet.");
    }
    @Override
    public void setReverse(boolean reverse) {
      throw new UnsupportedOperationException("Not supported yet.");
    }
    @Override
    public void templateCall(int endPosition, long hash, int index) {
      throw new UnsupportedOperationException("Not supported yet.");
    }
    @Override
    public TemplateCall threadClone(HashingRegion region) {
      throw new UnsupportedOperationException("Not supported yet.");
    }
    @Override
    public void threadFinish() {
      throw new UnsupportedOperationException("Not supported yet.");
    }
  }
  private static class FakeHashFunction implements NgsHashFunction {
    @Override
    public NgsHashFunction threadClone(HashingRegion region) {
      throw new UnsupportedOperationException("Not supported yet.");
    }
    @Override
    public void threadFinish() {
      throw new UnsupportedOperationException("Not supported yet.");
    }
    @Override
    public void readAll(int readId, boolean reverse) {
      throw new UnsupportedOperationException("Not supported yet.");
    }
    @Override
    public void setReadSequences(long numberReads) {
      throw new UnsupportedOperationException("Not supported yet.");
    }
    @Override
    public void setValues(int id2, boolean reverse) {
      throw new UnsupportedOperationException("Not supported yet.");
    }
    @Override
    public void hashStep(byte code) {
      throw new UnsupportedOperationException("Not supported yet.");
    }
    @Override
    public void hashStep() {
      throw new UnsupportedOperationException("Not supported yet.");
    }
    @Override
    public int numberWindows() {
      throw new UnsupportedOperationException("Not supported yet.");
    }
    @Override
    public int readLength() {
      throw new UnsupportedOperationException("Not supported yet.");
    }
    @Override
    public void reset() {
      throw new UnsupportedOperationException("Not supported yet.");
    }
    @Override
    public int windowSize() {
      throw new UnsupportedOperationException("Not supported yet.");
    }
    @Override
    public void endSequence() {
      throw new UnsupportedOperationException("Not supported yet.");
    }
    @Override
    public int fastScore(int readId) {
      throw new UnsupportedOperationException("Not supported yet.");
    }
    @Override
    public int indelScore(int readId) {
      throw new UnsupportedOperationException("Not supported yet.");
    }
    @Override
    public void logStatistics() {
      throw new UnsupportedOperationException("Not supported yet.");
    }
    @Override
    public void templateBidirectional(int endPosition) {
      throw new UnsupportedOperationException("Not supported yet.");
    }
    @Override
    public void templateForward(int endPosition) {
      throw new UnsupportedOperationException("Not supported yet.");
    }
    @Override
    public void templateReverse(int endPosition) {
      throw new UnsupportedOperationException("Not supported yet.");
    }
    @Override
    public void templateSet(long name, int length) {
      throw new UnsupportedOperationException("Not supported yet.");
    }
  }
}
