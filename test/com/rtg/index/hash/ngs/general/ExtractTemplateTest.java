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
package com.rtg.index.hash.ngs.general;


import com.rtg.index.hash.ngs.NgsHashFunction;
import com.rtg.index.hash.ngs.TemplateCall;
import com.rtg.launcher.HashingRegion;
import com.rtg.util.integrity.IntegralAbstract;

import junit.framework.TestCase;

/**
 */
public class ExtractTemplateTest extends TestCase {

  private static class TemplateCallMock extends IntegralAbstract implements TemplateCall, Cloneable {
    @Override
    public void templateCall(final int id, final long hash, final int index) {
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
    /**
     */
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

  public void test() {
    final StringBuilder sb = new StringBuilder();
    final SingleMask ms = new SingleMask(4, 8, 0, 1);
    final TemplateCall template = new TemplateCallMock();
    final ExtractTemplate r = new ExtractTemplate(ms, template, 1);
    r.toString(sb);
    final String s = sb.toString();
    //    System.out.println(s);
    assertTrue(s.contains(" index=1 end=-1"));
    assertTrue(s.contains("ExtractTemplate templateCall="));
  }

}
