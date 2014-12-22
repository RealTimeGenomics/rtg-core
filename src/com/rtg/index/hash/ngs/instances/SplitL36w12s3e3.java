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
 * w = 12 s = 3 e = 3
 * Window size 12
 * Guaranteed to find all triple substitutions.
 * Guaranteed to find all triple indels.
 * <br>
 * <table>
 * <tr> <th>window</th> <th> id </th> <th> indel </th></tr>
 * <tr> <td><code>ooooxx</code></td>  <td> 0 </td> <td> 1 </td></tr>
 * <tr> <td><code>oooxxo</code></td>  <td> 1 </td> <td> 1 </td></tr>
 * <tr> <td><code>ooxxoo</code></td>  <td> 2 </td> <td> 1 </td></tr>
 * <tr> <td><code>oxxooo</code></td>  <td> 3 </td> <td> 1 </td></tr>
 * <tr> <td><code>xxoooo</code></td>  <td> 4 </td> <td> 1 </td></tr>
 * <tr> <td><code>oooxox</code></td>  <td> 5 </td> <td> 3 </td></tr>
 * <tr> <td><code>ooxoox</code></td>  <td> 6 </td> <td> 3 </td></tr>
 * <tr> <td><code>xoxooo</code></td>  <td> 7 </td> <td> 3 </td></tr>
 * <tr> <td></td>                     <td> 8 </td> <td> 14 </td></tr>
 * </table>
 */
public class SplitL36w12s3e3 extends ImplementHashFunction {

  private static final int NUMBER_BITS = 24;

  private static final int NUMBER_WINDOWS = 8;

  /**
   * A factory for the outer class.
   */
  public static final HashFunctionFactory FACTORY = new HashFunctionFactory() {
    @Override
    public NgsHashFunction create(final ReadCall readCall, final TemplateCall templateCall) {
      return new SplitL36w12s3e3(readCall, templateCall);
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
      return 12;
    }
  };

  private static final long MASK_12 =    0xFFF;
  private static final long MASK_12_6 = MASK_12 << 6;
  private static final long MASK_12_12 = MASK_12 << 12;
  private static final long MASK_12_18 = MASK_12 << 18;
  private static final long MASK_12_24 = MASK_12 << 24;

  private static final long MASK_6 =    0x3F;
  private static final long MASK_6_6 = MASK_6 << 6;
  private static final long MASK_6_12 = MASK_6 << 12;
  private static final long MASK_6_12I = MASK_6 << 13;
  private static final long MASK_6_12D = MASK_6 << 11;
  private static final long MASK_6_18 = MASK_6 << 18;
  private static final long MASK_6_18I = MASK_6 << 19;
  private static final long MASK_6_18D = MASK_6 << 17;
  private static final long MASK_6_30 = MASK_6 << 30;
  private static final long MASK_6_30I = MASK_6 << 31;
  private static final long MASK_6_30D = MASK_6 << 29;

  /**
   * @param readCall used in subclasses to process results of read hits.
   * @param templateCall used in subclasses to process results of template hits.
   */
  public SplitL36w12s3e3(final ReadCall readCall, final TemplateCall templateCall) {
    super(36, 12, readCall, templateCall);
    setHashFunction();
  }

  @Override
  public void readAll(final int readId, long v0, long v1) {
    if (mSoFar < 36) {
      return;
    }
//    System.err.println("  v0=" + Utils.toBits(v0));
//    System.err.println("  v1=" + Utils.toBits(v1));
//    System.err.println();
    final long l0a = (v0 & MASK_12) << 12;
    final long l0b = v1 & MASK_12;
    final long l0 = l0a | l0b;
    mReadCall.readCall(readId, hash(l0), 0);

    final long l1a = (v0 & MASK_12_6) << 6;
    final long l1b = (v1 & MASK_12_6) >>> 6;
    final long l1 = l1a | l1b;
    mReadCall.readCall(readId, hash(l1), 1);

    final long l2a = v0 & MASK_12_12;
    final long l2b = (v1 & MASK_12_12) >>> 12;
    final long l2 = l2a | l2b;
    mReadCall.readCall(readId, hash(l2), 2);

    final long l3a = (v0 & MASK_12_18) >>> 6;
    final long l3b = (v1 & MASK_12_18) >>> 18;
    final long l3 = l3a | l3b;
    mReadCall.readCall(readId, hash(l3), 3);

    final long l4a = (v0 & MASK_12_24) >>> 12;
    final long l4b = (v1 & MASK_12_24) >>> 24;
    final long l4 = l4a | l4b;
    mReadCall.readCall(readId, hash(l4), 4);

    final long l5aa = v0 & MASK_6;
    final long l5ab = (v0 & MASK_6_12) >>> 6;
    final long l5a = l5aa | l5ab;
    final long l5ba = (v1 & MASK_6) << 12;
    final long l5bb = (v1 & MASK_6_12) << 6;
    final long l5b = l5ba | l5bb;
    final long l5 = l5a | l5b;
    mReadCall.readCall(readId, hash(l5), 5);

    final long l6aa = (v0 & MASK_6_6) >>> 6;
    final long l6ab = (v0 & MASK_6_18) >>> 12;
    final long l6a = l6aa | l6ab;
    final long l6ba = (v1 & MASK_6_6) << 6;
    final long l6bb =  v1 & MASK_6_18;
    final long l6b = l6ba | l6bb;
    final long l6 = l6a | l6b;
    mReadCall.readCall(readId, hash(l6), 6);

    final long l7aa = (v0 & MASK_6_18) >>> 18;
    final long l7ab = (v0 & MASK_6_30) >>> 24;
    final long l7a = l7aa | l7ab;
    final long l7ba = (v1 & MASK_6_18) >>> 6;
    final long l7bb = (v1 & MASK_6_30) >>> 12;
    final long l7b = l7ba | l7bb;
    final long l7 = l7a | l7b;
    mReadCall.readCall(readId, hash(l7), 7);
}

  @Override
  public void templateAll(final int endPosition, long v0, long v1) throws IOException {
    if (mSoFar < 36) {
      return;
    }
    //System.err.println("  v0=" + Utils.toBitsSep(v0));
    //System.err.println("  v1=" + Utils.toBitsSep(v1));
    //System.err.println();

    final long l0a = (v0 & MASK_12) << 12;
    final long l0b =  v1 & MASK_12;
    final long l0 = l0a | l0b;
    mTemplateCall.templateCall(endPosition, hash(l0), 0);

    final long l1a = (v0 & MASK_12_6) << 6;
    final long l1b = (v1 & MASK_12_6) >>> 6;
    final long l1 = l1a | l1b;
    mTemplateCall.templateCall(endPosition, hash(l1), 1);

    final long l2a = v0 & MASK_12_12;
    final long l2b = (v1 & MASK_12_12) >>> 12;
    final long l2 = l2a | l2b;
    mTemplateCall.templateCall(endPosition, hash(l2), 2);

    final long l3a = (v0 & MASK_12_18) >>> 6;
    final long l3b = (v1 & MASK_12_18) >>> 18;
    final long l3 = l3a | l3b;
    mTemplateCall.templateCall(endPosition, hash(l3), 3);

    final long l4a = (v0 & MASK_12_24) >>> 12;
    final long l4b = (v1 & MASK_12_24) >>> 24;
    final long l4 = l4a | l4b;
    mTemplateCall.templateCall(endPosition, hash(l4), 4);

    final long l5aa = v0 & MASK_6;
    final long l5ab = (v0 & MASK_6_12) >>> 6;
    final long l5a = l5aa | l5ab;
    final long l5ba = (v1 & MASK_6) << 12;
    final long l5bb = (v1 & MASK_6_12) << 6;
    final long l5b = l5ba | l5bb;
    final long l5 = l5a | l5b;
    mTemplateCall.templateCall(endPosition, hash(l5), 5);

    //insert
    final long l5abi = (v0 & MASK_6_12I) >>> 7;
    final long l5ai = l5aa | l5abi;
    final long l5bbi = (v1 & MASK_6_12I) << 5;
    final long l5bi = l5ba | l5bbi;
    final long l5i = l5ai | l5bi;
    mTemplateCall.templateCall(endPosition, hash(l5i), 5);

    //delete
    final long l5abd = (v0 & MASK_6_12D) >>> 5;
    final long l5ad = l5aa | l5abd;
    final long l5bbd = (v1 & MASK_6_12D) << 7;
    final long l5bd = l5ba | l5bbd;
    final long l5d = l5ad | l5bd;
    mTemplateCall.templateCall(endPosition, hash(l5d), 5);

    final long l6aa = (v0 & MASK_6_6) >>> 6;
    final long l6ab = (v0 & MASK_6_18) >>> 12;
    final long l6a = l6aa | l6ab;
    final long l6ba = (v1 & MASK_6_6) << 6;
    final long l6bb =  v1 & MASK_6_18;
    final long l6b = l6ba | l6bb;
    final long l6 = l6a | l6b;
    mTemplateCall.templateCall(endPosition, hash(l6), 6);

    //insert
    final long l6abi = (v0 & MASK_6_18I) >>> 13;
    final long l6ai = l6aa | l6abi;
    final long l6bbi = (v1 & MASK_6_18I) >>> 1;
    final long l6bi = l6ba | l6bbi;
    final long l6i = l6ai | l6bi;
    mTemplateCall.templateCall(endPosition, hash(l6i), 6);

    //delete
    final long l6abd = (v0 & MASK_6_18D) >>> 11;
    final long l6ad = l6aa | l6abd;
    final long l6bbd = (v1 & MASK_6_18D) << 1;
    final long l6bd = l6ba | l6bbd;
    final long l6d = l6ad | l6bd;
    mTemplateCall.templateCall(endPosition, hash(l6d), 6);

    final long l7aa = (v0 & MASK_6_18) >>> 18;
    final long l7ab = (v0 & MASK_6_30) >>> 24;
    final long l7a = l7aa | l7ab;
    final long l7ba = (v1 & MASK_6_18) >>> 6;
    final long l7bb = (v1 & MASK_6_30) >>> 12;
    final long l7b = l7ba | l7bb;
    final long l7 = l7a | l7b;
    mTemplateCall.templateCall(endPosition, hash(l7), 7);

    //insert
    final long l7abi = (v0 & MASK_6_30I) >>> 25;
    final long l7ai = l7aa | l7abi;
    final long l7bbi = (v1 & MASK_6_30I) >>> 13;
    final long l7bi = l7ba | l7bbi;
    final long l7i = l7ai | l7bi;
    mTemplateCall.templateCall(endPosition, hash(l7i), 7);

    //delete
    final long l7abd = (v0 & MASK_6_30D) >>> 23;
    final long l7ad = l7aa | l7abd;
    final long l7bbd = (v1 & MASK_6_30D) >>> 11;
    final long l7bd = l7ba | l7bbd;
    final long l7d = l7ad | l7bd;
    mTemplateCall.templateCall(endPosition, hash(l7d), 7);

    mTemplateCall.done();
  }

  @Override
  public int numberWindows() {
    return NUMBER_WINDOWS;
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("Split l=").append(readLength()).append(" w=").append(windowSize()).append(" e=3");
  }
}

