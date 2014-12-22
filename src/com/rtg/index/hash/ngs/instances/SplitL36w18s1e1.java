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
 * w=18 s=1 e=1
 * Guaranteed to find all single substitutions, single indels.
 *
 * Two split windows with the following patterns (blocks of size 18).
 * Does not explicitly handle indels.
 * <br>
 * <table>
 * <tr> <th>window</th> <th> id </th> <th> indel </th> </tr>
 * <tr> <td><code>ox</code></td>  <td> 0 </td> <td> 1 </td></tr>
 * <tr> <td><code>xo</code></td>  <td> 1 </td> <td> 1 </td></tr>
 * <tr> <td></td>                 <td> 2 </td> <td> 2 </td></tr>
 * </table>
 */
public class SplitL36w18s1e1 extends ImplementHashFunction {

  private static final int NUMBER_BITS = 36;

  private static final int NUMBER_WINDOWS = 2;

  /**
   * A factory for the outer class.
   */
  public static final HashFunctionFactory FACTORY = new HashFunctionFactory() {
    @Override
    public NgsHashFunction create(final ReadCall readCall, final TemplateCall templateCall) {
      return new SplitL36w18s1e1(readCall, templateCall);
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
      return 18;
    }
  };

  private static final long MASK_18 =    0x3FFFF;
  private static final long MASK_18_18 = MASK_18 << 18;

  /**
   * @param readCall used in subclasses to process results of read hits.
   * @param templateCall used in subclasses to process results of template hits.
   */
  public SplitL36w18s1e1(final ReadCall readCall, final TemplateCall templateCall) {
    super(36, 18, readCall, templateCall);
    setHashFunction();
  }

  @Override
  public void readAll(final int readId, long v0, long v1) {
    if (mSoFar < 36) {
      return;
    }
    final long la = (v0 & MASK_18) << 18;
    final long lb =  v1 & MASK_18;
    final long l2x0 = la | lb;
    mReadCall.readCall(readId, hash(l2x0), 0);

    final long lm = v0 & MASK_18_18;
    final long ln = (v1 & MASK_18_18) >>> 18;
    final long l2x2 = lm | ln;
    mReadCall.readCall(readId, hash(l2x2), 1);
  }

  @Override
  public void templateAll(final int endPosition, long v0, long v1) throws IOException {
    if (mSoFar < 36) {
      return;
    }
    final long la = (v0 & MASK_18) << 18;
    final long lb =  v1 & MASK_18;
    final long l2x0 = la | lb;
    mTemplateCall.templateCall(endPosition, hash(l2x0), 0);

    final long lm = v0 & MASK_18_18;
    final long ln = (v1 & MASK_18_18) >>> 18;
    final long l2x2 = lm | ln;
    mTemplateCall.templateCall(endPosition, hash(l2x2), 1);

    mTemplateCall.done();
  }

  @Override
  public int numberWindows() {
    return NUMBER_WINDOWS;
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("Split l=").append(readLength()).append(" w=").append(windowSize()).append(" s=1");
  }
}

