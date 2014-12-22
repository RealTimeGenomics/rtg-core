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
 * Automatically generated from masks file.
 * DO NOT MODIFY.
 *
 * l = 51 w = 18 s = 2 e = 1
 *
 * <br>
 * <table border="1">
 * <tr> <th>window</th> <th> id </th> <th> indel </th></tr>
 * <tr>  <td> 00000111</td> <td> 0</td> <td> 1</td> </tr>
 * <tr>  <td> 00011100</td> <td> 1</td> <td> 1</td> </tr>
 * <tr>  <td> 00111000</td> <td> 2</td> <td> 1</td> </tr>
 * <tr>  <td> 11100000</td> <td> 3</td> <td> 1</td> </tr>
 * <tr>  <td> 00001011</td> <td> 4</td> <td> 1</td> </tr>
 * <tr>  <td> </td> <td> 5</td> <td> 5</td> </tr>
 * </table>
 *
 */
public class MaskL51w18s2e1 extends ImplementHashFunction {

  private static final int LENGTH = 51;

  private static final int NUMBER_NT = 18;

  private static final int NUMBER_BITS = 2 * NUMBER_NT;

  private static final int NUMBER_WINDOWS = 5;

  /**
   * A factory for the outer class.
   */
  public static final HashFunctionFactory FACTORY = new HashFunctionFactory() {
    @Override
    public NgsHashFunction create(final ReadCall readCall, final TemplateCall templateCall) {
      return new MaskL51w18s2e1(readCall, templateCall);
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
      return NUMBER_NT;
    }
  };

  /**
   * @param readCall used in subclasses to process results of read hits.
   * @param templateCall used in subclasses to process results of template hits.
   */
  public MaskL51w18s2e1(final ReadCall readCall, final TemplateCall templateCall) {
    super(LENGTH, NUMBER_NT, readCall, templateCall);
    setHashFunction();
  }

  @Override
  public void readAll(final int readId, long v0, long v1) {
    if (mSoFar < LENGTH) {
      return;
    }

    final long l0 =
    (v0 & ((1L << 18) - 1L)) | ((v1 & ((1L << 18) - 1L)) << 18);
    mReadCall.readCall(readId, hash(l0), 0);

    final long l1 =
    ((v0 >> 12) & ((1L << 18) - 1L)) | (((v1 >> 12) & ((1L << 18) - 1L)) << 18);
    mReadCall.readCall(readId, hash(l1), 1);

    final long l2 =
    ((v0 >> 18) & ((1L << 18) - 1L)) | (((v1 >> 18) & ((1L << 18) - 1L)) << 18);
    mReadCall.readCall(readId, hash(l2), 2);

    final long l3 =
    ((v0 >> 30) & ((1L << 18) - 1L)) | (((v1 >> 30) & ((1L << 18) - 1L)) << 18);
    mReadCall.readCall(readId, hash(l3), 3);

    final long l4 =
    (v0 & ((1L << 12) - 1L)) | ((v1 & ((1L << 12) - 1L)) << 18)
    | (((v0 >> 18) & ((1L << 6) - 1L)) << 12) | (((v1 >> 18) & ((1L << 6) - 1L)) << 30);
    mReadCall.readCall(readId, hash(l4), 4);
}

  @Override
  public void templateAll(final int endPosition, long v0, long v1) throws IOException {
    if (mSoFar < LENGTH) {
      return;
    }

    final long l0 =
    (v0 & ((1L << 18) - 1L)) | ((v1 & ((1L << 18) - 1L)) << 18);
    mTemplateCall.templateCall(endPosition, hash(l0), 0);

    final long l1 =
    ((v0 >> 12) & ((1L << 18) - 1L)) | (((v1 >> 12) & ((1L << 18) - 1L)) << 18);
    mTemplateCall.templateCall(endPosition, hash(l1), 1);

    final long l2 =
    ((v0 >> 18) & ((1L << 18) - 1L)) | (((v1 >> 18) & ((1L << 18) - 1L)) << 18);
    mTemplateCall.templateCall(endPosition, hash(l2), 2);

    final long l3 =
    ((v0 >> 30) & ((1L << 18) - 1L)) | (((v1 >> 30) & ((1L << 18) - 1L)) << 18);
    mTemplateCall.templateCall(endPosition, hash(l3), 3);

    final long l4 =
    (v0 & ((1L << 12) - 1L)) | ((v1 & ((1L << 12) - 1L)) << 18)
    | (((v0 >> 18) & ((1L << 6) - 1L)) << 12) | (((v1 >> 18) & ((1L << 6) - 1L)) << 30);
    mTemplateCall.templateCall(endPosition, hash(l4), 4);

    mTemplateCall.done();
  }

  @Override
  public int numberWindows() {
    return NUMBER_WINDOWS;
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("MaskL51w18s2e1 l=" + LENGTH + " w=" + NUMBER_NT + " s=2 e=1");
  }
}
