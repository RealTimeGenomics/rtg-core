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
 * w=18 s=2 e=1
 * Guaranteed to find all two substitutions, single indels.
 *
 * Six split windows with the following patterns (blocks of size 9).
 * Does not handle indels.
 * <br>
 * <table summary="mask cases">
 * <tr> <th>window</th> <th> id </th> <th> indel </th></tr>
 * <tr> <td><code>ooxx</code></td>  <td> 0 </td> <td> 1 </td></tr>
 * <tr> <td><code>oxxo</code></td>  <td> 1 </td> <td> 1 </td></tr>
 * <tr> <td><code>xxoo</code></td>  <td> 2 </td> <td> 1 </td></tr>
 * <tr> <td><code>oxox</code></td>  <td> 3 </td> <td> 1 </td></tr>
 * <tr> <td><code>xoxo</code></td>  <td> 4 </td> <td> 1 </td></tr>
 * <tr> <td><code>xoox</code></td>  <td> 5 </td> <td> 1 </td></tr>
 * <tr> <td></td>                   <td> 6 </td> <td> 6 </td></tr>
 * </table>
 */
public class SplitL36w18s2e1 extends ImplementHashFunction {

  private static final int NUMBER_BITS = 36;

  private static final int NUMBER_WINDOWS = 6;

  /**
   * A factory for the outer class.
   */
  public static final HashFunctionFactory FACTORY = new HashFunctionFactory() {
    @Override
    public NgsHashFunction create(final ReadCall readCall, final TemplateCall templateCall) {
      return new SplitL36w18s2e1(readCall, templateCall);
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

  private static final long MASK_9 =     0x1FF;
  private static final long MASK_9_27 =  MASK_9 << 27;
  private static final long MASK_9_18 =  MASK_9 << 18;
  private static final long MASK_9_9 =   MASK_9 << 9;
  private static final long MASK_18 =    0x3FFFF;
  private static final long MASK_18_18 = MASK_18 << 18;
  private static final long MASK_18_9 =  MASK_18 << 9;

  /**
   * @param readCall used in subclasses to process results of read hits.
   * @param templateCall used in subclasses to process results of template hits.
   */
  public SplitL36w18s2e1(final ReadCall readCall, final TemplateCall templateCall) {
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

    final long lc = (v0 & MASK_18_9) << 9;
    final long ld = (v1 >>> 9) & MASK_18;
    final long l2x1 = lc | ld;
    mReadCall.readCall(readId, hash(l2x1), 1);

    final long lm = v0 & MASK_18_18;
    final long ln = (v1 & MASK_18_18) >>> 18;
    final long l2x2 = lm | ln;
    mReadCall.readCall(readId, hash(l2x2), 2);

    final long le = (v0 & MASK_9_18) << 9;
    final long lf = v1 & MASK_9_18;
    final long lg = (v0 & MASK_9) << 9;
    final long lh = v1 & MASK_9;
    final long l1x0 = le | lf | lg | lh;
    mReadCall.readCall(readId, hash(l1x0), 3);

    final long li = v0 & MASK_9_27;
    final long lj = (v1 & MASK_9_27) >>> 9;
    final long lk = v0 & MASK_9_9;
    final long ll = (v1 & MASK_9_9) >>> 9;
    final long l1x1 = li | lj | lk | ll;
    mReadCall.readCall(readId, hash(l1x1), 4);

    final long l0x0 = li | lj | lg | lh;
    mReadCall.readCall(readId, hash(l0x0), 5);
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

    final long lc = (v0 & MASK_18_9) << 9;
    final long ld = (v1 >>> 9) & MASK_18;
    final long l2x1 = lc | ld;
    mTemplateCall.templateCall(endPosition, hash(l2x1), 1);

    final long lm = v0 & MASK_18_18;
    final long ln = (v1 & MASK_18_18) >>> 18;
    final long l2x2 = lm | ln;
    mTemplateCall.templateCall(endPosition, hash(l2x2), 2);

    final long le = (v0 & MASK_9_18) << 9;
    final long lf = v1 & MASK_9_18;
    final long lg = (v0 & MASK_9) << 9;
    final long lh = v1 & MASK_9;
    final long l1x0 = le | lf | lg | lh;
    mTemplateCall.templateCall(endPosition, hash(l1x0), 3);

    final long li = v0 & MASK_9_27;
    final long lj = (v1 & MASK_9_27) >>> 9;
    final long lk = v0 & MASK_9_9;
    final long ll = (v1 & MASK_9_9) >>> 9;
    final long l1x1 = li | lj | lk | ll;
    mTemplateCall.templateCall(endPosition, hash(l1x1), 4);

    final long l0x0 = li | lj | lg | lh;
    mTemplateCall.templateCall(endPosition, hash(l0x0), 5);

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

