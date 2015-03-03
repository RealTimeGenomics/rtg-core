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
import com.rtg.index.hash.ngs.ImplementHashFunction;
import com.rtg.index.hash.ngs.NgsHashFunction;
import com.rtg.index.hash.ngs.ReadCall;
import com.rtg.index.hash.ngs.TemplateCall;

/**
 * w=4 s=0 e=0
 * Single split windows with the following patterns.
 * Does not handle indels.
 * Intended for testing.
 * <br>
 * <table summary="mask cases">
 * <tr> <th>window</th> <th> id </th> </tr>
 * <tr> <td><code>xxxx</code></td>  <td> 0 </td> </tr>
 * </table>
 */
public class SplitL4w4s0e0 extends ImplementHashFunction {

  private static final int NUMBER_BITS = 8;

  private static final int NUMBER_WINDOWS = 1;

  private static final int MASK_4 = 0xF;

  /**
   * A factory for the outer class.
   */
  public static final HashFunctionFactory FACTORY = new HashFunctionFactory() {
    @Override
    public NgsHashFunction create(final ReadCall readCall, final TemplateCall templateCall) {
      return new SplitL4w4s0e0(readCall, templateCall);
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
      return 4;
    }
  };

  /**
   * @param readCall used in subclasses to process results of read hits.
   * @param templateCall used in subclasses to process results of template hits.
   */
  public SplitL4w4s0e0(final ReadCall readCall, final TemplateCall templateCall) {
    super(4, 4, readCall, templateCall);
    setHashFunction();
  }

  @Override
  public void readAll(final int readId, long v0, long v1) {
    if (mSoFar < 4) {
      return;
    }

    final long la = (v0 & MASK_4) << 4;
    final long lb = v1 & MASK_4;
    final long l = la | lb;
    mReadCall.readCall(readId, hash(l), 0);
  }

  @Override
  public void templateAll(final int endPosition, long v0, long v1) throws IOException {
    if (mSoFar < 4) {
      return;
    }

    final long la = (v0 & MASK_4) << 4;
    final long lb = v1 & MASK_4;
    final long l = la | lb;
    mTemplateCall.templateCall(endPosition, hash(l), 0);

    mTemplateCall.done();
  }

  @Override
  public int numberWindows() {
    return NUMBER_WINDOWS;
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("Split l=").append(readLength()).append(" w=").append(windowSize()).append(" no indels");
  }
}

