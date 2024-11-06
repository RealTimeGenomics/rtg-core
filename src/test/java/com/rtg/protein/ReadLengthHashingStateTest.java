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
