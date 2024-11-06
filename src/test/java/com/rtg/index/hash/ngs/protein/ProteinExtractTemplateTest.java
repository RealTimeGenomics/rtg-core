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
package com.rtg.index.hash.ngs.protein;

import java.io.IOException;

import com.rtg.index.hash.ngs.NgsHashFunction;
import com.rtg.index.hash.ngs.TemplateCall;
import com.rtg.index.hash.ngs.general.SingleMask;
import com.rtg.launcher.HashingRegion;
import com.rtg.util.integrity.IntegralAbstract;

import junit.framework.TestCase;
public class ProteinExtractTemplateTest extends TestCase {

  private static class TemplateCallMock extends IntegralAbstract implements TemplateCall, Cloneable {

    @Override
    public void templateCall(final int id, final long hash, final int index) {
      if (hash == 0) {
        fail();
      }
    }

    @Override
    public void done() {
    }

    @Override
    public void endSequence() {
    }

    @Override
    public void set(final long name, final int length) {
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
    }
    @Override
    public TemplateCall clone() throws CloneNotSupportedException {
      return (TemplateCall) super.clone();
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
    public void toString(final StringBuilder sb) {
    }

    @Override
    public boolean integrity() {
      return true;
    }

    @Override
    public void logStatistics() {
      // do nothing
    }
  }

  public void test() throws IOException {
    final StringBuilder sb = new StringBuilder();
    final SingleMask ms = new SingleMask(4, 8, 0, 1);
    final TemplateCall template = new TemplateCallMock();
    final ProteinExtractTemplate r = new ProteinExtractTemplate(ms, template, 1);
    assertTrue(r.integrity());
    r.masked(0);
    r.toString(sb);
    final String s = sb.toString();
    assertTrue(s.contains(" index=1 end=-1"));
    assertTrue(s.contains("ProteinExtractTemplate templateCall="));
  }
}
