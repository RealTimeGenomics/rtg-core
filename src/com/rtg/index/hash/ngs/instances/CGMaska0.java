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
package com.rtg.index.hash.ngs.instances;

import java.io.IOException;

import com.rtg.index.hash.ngs.HashFunctionFactory;
import com.rtg.index.hash.ngs.NgsHashFunction;
import com.rtg.index.hash.ngs.ReadCall;
import com.rtg.index.hash.ngs.TemplateCall;

/**
 * r=35 s=0
 * Doesn't handle much - here mainly for testing CG stuff before
 * implementing some real masks.
 */
public class CGMaska0 extends AbstractCGMask {

  private static final int READ_LENGTH = 35;

  private static final int WINDOW_LENGTH = 18;

  private static final int NUMBER_BITS = 36;

  private static final int NUMBER_WINDOWS = 2;

  /**
   * A factory for the outer class.
   */
  public static final HashFunctionFactory FACTORY = new CGHashFunctionFactory() {
    @Override
    public NgsHashFunction create(final ReadCall readCall, final TemplateCall templateCall) {
      return new CGMaska0(readCall, templateCall);
    }

    @Override
    public int hashBits() {
      return NUMBER_BITS;
    }

    @Override
    public int windowBits() {
      return NUMBER_BITS;
    }

    @Override
    public int numberWindows() {
      return NUMBER_WINDOWS;
    }

    @Override
    public int windowSize() {
      return WINDOW_LENGTH;
    }
  };

  private static final long MASK_A = ((1L << 18) - 1) << 17;
  private static final long MASK_B = (1L << 17) - 1;

  /**
   * @param readCall used in subclasses to process results of read hits.
   * @param templateCall used in subclasses to process results of template hits.
   */
  public CGMaska0(final ReadCall readCall, final TemplateCall templateCall) {
    super(READ_LENGTH, WINDOW_LENGTH, readCall, templateCall, LONG_BITS - (READ_LENGTH + 6));
    setHashFunction();
  }

  @Override
  public void readAll(final int readId, long v0, long v1) {
    //System.err.println("readAll sofar=" + mSoFar);
    if (mSoFar < READ_LENGTH) {
      return;
    }
    final long la = (v0 & MASK_A) << 1;
    final long lb = (v1 & MASK_A) >>> 17;
    final long l2x0 = la | lb;
    mReadCall.readCall(readId, hash(l2x0), 0);

    final long lc = (v0 & MASK_B) << 18;
    final long ld =  v1 & MASK_B;
    final long l2x1 = lc | ld;
    mReadCall.readCall(readId, hash(l2x1), 1);
  }

  private static final long MASKT_A = ((1L << 18) - 1) << (17 + 6);
  private static final long MASKT_B0 = (1L << 10) - 1;
  private static final long MASKT_B1 = ((1L << 7) - 1) << 16;

  @Override
  public void templateAll(final int endPosition, long v0, long v1) throws IOException {
    //System.err.println("templateAll sofar=" + mSoFar + " endPosition=" + endPosition);
    //System.err.println("    v0=" + Utils.toBitsSep(v0));
    //System.err.println("    v1=" + Utils.toBitsSep(v1));
    if (mSoFar < READ_LENGTH + 1) {
      return;
    }

    //The -7 allows for an insert of 6nt the +2 is experimenntal
    final int endP = endPosition - 6 + 2;

    //System.err.println("templateAll index 0");
    final long la = (v0 & MASKT_A) >>> 5;
    final long lb = (v1 & MASKT_A) >>> 23;
    final long l2x0 = la | lb;
    mTemplateCall.templateCall(endP, hash(l2x0), 0);

    //System.err.println("templateAll index 1");
    final long lc0 = (v0 & MASKT_B0) << 18;
    final long lc1 = (v0 & MASKT_B1) << 12;
    final long lc = lc0 | lc1;
    final long ld0 =  v1 & MASKT_B0;
    final long ld1 = (v1 & MASKT_B1) >>> 6;
    final long ld = ld0 | ld1;
    final long l2x1 = lc | ld;
    mTemplateCall.templateCall(endP, hash(l2x1), 1);

    mTemplateCall.done();
  }

  @Override
  public int numberWindows() {
    return NUMBER_WINDOWS;
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("CG l=").append(readLength()).append(" w=").append(windowSize()).append(" e=0");
  }
}

