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
 * l = 36 w = 18 s = 3 e = 1
 *
 * <br>
 * <table border="1">
 * <tr> <th>window</th> <th> id </th> <th> indel </th></tr>
 * <tr>  <td> 000111</td> <td> 0</td> <td> 1</td> </tr>
 * <tr>  <td> 001110</td> <td> 1</td> <td> 1</td> </tr>
 * <tr>  <td> 011100</td> <td> 2</td> <td> 1</td> </tr>
 * <tr>  <td> 111000</td> <td> 3</td> <td> 1</td> </tr>
 * <tr>  <td> 001011</td> <td> 4</td> <td> 1</td> </tr>
 * <tr>  <td> 010110</td> <td> 5</td> <td> 1</td> </tr>
 * <tr>  <td> 101100</td> <td> 6</td> <td> 1</td> </tr>
 * <tr>  <td> 010011</td> <td> 7</td> <td> 1</td> </tr>
 * <tr>  <td> 100110</td> <td> 8</td> <td> 1</td> </tr>
 * <tr>  <td> 100011</td> <td> 9</td> <td> 1</td> </tr>
 * <tr>  <td> 001101</td> <td> 10</td> <td> 1</td> </tr>
 * <tr>  <td> 011010</td> <td> 11</td> <td> 1</td> </tr>
 * <tr>  <td> 110100</td> <td> 12</td> <td> 1</td> </tr>
 * <tr>  <td> 010101</td> <td> 13</td> <td> 1</td> </tr>
 * <tr>  <td> 101010</td> <td> 14</td> <td> 1</td> </tr>
 * <tr>  <td> 100101</td> <td> 15</td> <td> 1</td> </tr>
 * <tr>  <td> 011001</td> <td> 16</td> <td> 1</td> </tr>
 * <tr>  <td> 110010</td> <td> 17</td> <td> 1</td> </tr>
 * <tr>  <td> 101001</td> <td> 18</td> <td> 1</td> </tr>
 * <tr>  <td> 110001</td> <td> 19</td> <td> 1</td> </tr>
 * <tr>  <td> </td> <td> 20</td> <td> 20</td> </tr>
 * </table>
 *
 */
public class MaskL36w18s3e1 extends ImplementHashFunction {

  private static final int LENGTH = 36;

  private static final int NUMBER_NT = 18;

  private static final int NUMBER_BITS = 2 * NUMBER_NT;

  private static final int NUMBER_WINDOWS = 20;

  /**
   * A factory for the outer class.
   */
  public static final HashFunctionFactory FACTORY = new HashFunctionFactory() {
    @Override
    public NgsHashFunction create(final ReadCall readCall, final TemplateCall templateCall) {
      return new MaskL36w18s3e1(readCall, templateCall);
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
  public MaskL36w18s3e1(final ReadCall readCall, final TemplateCall templateCall) {
    super(LENGTH, NUMBER_NT, readCall, templateCall);
    setHashFunction();
  }

  @Override
  public void readAll(final int readId, long v0, long v1) {
    final long l0 =
    (v0 & ((1L << 18) - 1L)) | ((v1 & ((1L << 18) - 1L)) << 18);
    mReadCall.readCall(readId, hash(l0), 0);

    final long l1 =
    ((v0 >> 6) & ((1L << 18) - 1L)) | (((v1 >> 6) & ((1L << 18) - 1L)) << 18);
    mReadCall.readCall(readId, hash(l1), 1);

    final long l2 =
    ((v0 >> 12) & ((1L << 18) - 1L)) | (((v1 >> 12) & ((1L << 18) - 1L)) << 18);
    mReadCall.readCall(readId, hash(l2), 2);

    final long l3 =
    ((v0 >> 18) & ((1L << 18) - 1L)) | (((v1 >> 18) & ((1L << 18) - 1L)) << 18);
    mReadCall.readCall(readId, hash(l3), 3);

    final long l4 =
    (v0 & ((1L << 12) - 1L)) | ((v1 & ((1L << 12) - 1L)) << 18)
    | (((v0 >> 18) & ((1L << 6) - 1L)) << 12) | (((v1 >> 18) & ((1L << 6) - 1L)) << 30);
    mReadCall.readCall(readId, hash(l4), 4);

    final long l5 =
    ((v0 >> 6) & ((1L << 12) - 1L)) | (((v1 >> 6) & ((1L << 12) - 1L)) << 18)
    | (((v0 >> 24) & ((1L << 6) - 1L)) << 12) | (((v1 >> 24) & ((1L << 6) - 1L)) << 30);
    mReadCall.readCall(readId, hash(l5), 5);

    final long l6 =
    ((v0 >> 12) & ((1L << 12) - 1L)) | (((v1 >> 12) & ((1L << 12) - 1L)) << 18)
    | (((v0 >> 30) & ((1L << 6) - 1L)) << 12) | (((v1 >> 30) & ((1L << 6) - 1L)) << 30);
    mReadCall.readCall(readId, hash(l6), 6);

    final long l7 =
    (v0 & ((1L << 12) - 1L)) | ((v1 & ((1L << 12) - 1L)) << 18)
    | (((v0 >> 24) & ((1L << 6) - 1L)) << 12) | (((v1 >> 24) & ((1L << 6) - 1L)) << 30);
    mReadCall.readCall(readId, hash(l7), 7);

    final long l8 =
    ((v0 >> 6) & ((1L << 12) - 1L)) | (((v1 >> 6) & ((1L << 12) - 1L)) << 18)
    | (((v0 >> 30) & ((1L << 6) - 1L)) << 12) | (((v1 >> 30) & ((1L << 6) - 1L)) << 30);
    mReadCall.readCall(readId, hash(l8), 8);

    final long l9 =
    (v0 & ((1L << 12) - 1L)) | ((v1 & ((1L << 12) - 1L)) << 18)
    | (((v0 >> 30) & ((1L << 6) - 1L)) << 12) | (((v1 >> 30) & ((1L << 6) - 1L)) << 30);
    mReadCall.readCall(readId, hash(l9), 9);

    final long l10 =
    (v0 & ((1L << 6) - 1L)) | ((v1 & ((1L << 6) - 1L)) << 18)
    | (((v0 >> 12) & ((1L << 12) - 1L)) << 6) | (((v1 >> 12) & ((1L << 12) - 1L)) << 24);
    mReadCall.readCall(readId, hash(l10), 10);

    final long l11 =
    ((v0 >> 6) & ((1L << 6) - 1L)) | (((v1 >> 6) & ((1L << 6) - 1L)) << 18)
    | (((v0 >> 18) & ((1L << 12) - 1L)) << 6) | (((v1 >> 18) & ((1L << 12) - 1L)) << 24);
    mReadCall.readCall(readId, hash(l11), 11);

    final long l12 =
    ((v0 >> 12) & ((1L << 6) - 1L)) | (((v1 >> 12) & ((1L << 6) - 1L)) << 18)
    | (((v0 >> 24) & ((1L << 12) - 1L)) << 6) | (((v1 >> 24) & ((1L << 12) - 1L)) << 24);
    mReadCall.readCall(readId, hash(l12), 12);

    final long l13 =
    (v0 & ((1L << 6) - 1L)) | ((v1 & ((1L << 6) - 1L)) << 18)
    | (((v0 >> 12) & ((1L << 6) - 1L)) << 6) | (((v1 >> 12) & ((1L << 6) - 1L)) << 24)
    | (((v0 >> 24) & ((1L << 6) - 1L)) << 12) | (((v1 >> 24) & ((1L << 6) - 1L)) << 30);
    mReadCall.readCall(readId, hash(l13), 13);

    final long l14 =
    ((v0 >> 6) & ((1L << 6) - 1L)) | (((v1 >> 6) & ((1L << 6) - 1L)) << 18)
    | (((v0 >> 18) & ((1L << 6) - 1L)) << 6) | (((v1 >> 18) & ((1L << 6) - 1L)) << 24)
    | (((v0 >> 30) & ((1L << 6) - 1L)) << 12) | (((v1 >> 30) & ((1L << 6) - 1L)) << 30);
    mReadCall.readCall(readId, hash(l14), 14);

    final long l15 =
    (v0 & ((1L << 6) - 1L)) | ((v1 & ((1L << 6) - 1L)) << 18)
    | (((v0 >> 12) & ((1L << 6) - 1L)) << 6) | (((v1 >> 12) & ((1L << 6) - 1L)) << 24)
    | (((v0 >> 30) & ((1L << 6) - 1L)) << 12) | (((v1 >> 30) & ((1L << 6) - 1L)) << 30);
    mReadCall.readCall(readId, hash(l15), 15);

    final long l16 =
    (v0 & ((1L << 6) - 1L)) | ((v1 & ((1L << 6) - 1L)) << 18)
    | (((v0 >> 18) & ((1L << 12) - 1L)) << 6) | (((v1 >> 18) & ((1L << 12) - 1L)) << 24);
    mReadCall.readCall(readId, hash(l16), 16);

    final long l17 =
    ((v0 >> 6) & ((1L << 6) - 1L)) | (((v1 >> 6) & ((1L << 6) - 1L)) << 18)
    | (((v0 >> 24) & ((1L << 12) - 1L)) << 6) | (((v1 >> 24) & ((1L << 12) - 1L)) << 24);
    mReadCall.readCall(readId, hash(l17), 17);

    final long l18 =
    (v0 & ((1L << 6) - 1L)) | ((v1 & ((1L << 6) - 1L)) << 18)
    | (((v0 >> 18) & ((1L << 6) - 1L)) << 6) | (((v1 >> 18) & ((1L << 6) - 1L)) << 24)
    | (((v0 >> 30) & ((1L << 6) - 1L)) << 12) | (((v1 >> 30) & ((1L << 6) - 1L)) << 30);
    mReadCall.readCall(readId, hash(l18), 18);

    final long l19 =
    (v0 & ((1L << 6) - 1L)) | ((v1 & ((1L << 6) - 1L)) << 18)
    | (((v0 >> 24) & ((1L << 12) - 1L)) << 6) | (((v1 >> 24) & ((1L << 12) - 1L)) << 24);
    mReadCall.readCall(readId, hash(l19), 19);
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
    ((v0 >> 6) & ((1L << 18) - 1L)) | (((v1 >> 6) & ((1L << 18) - 1L)) << 18);
    mTemplateCall.templateCall(endPosition, hash(l1), 1);

    final long l2 =
    ((v0 >> 12) & ((1L << 18) - 1L)) | (((v1 >> 12) & ((1L << 18) - 1L)) << 18);
    mTemplateCall.templateCall(endPosition, hash(l2), 2);

    final long l3 =
    ((v0 >> 18) & ((1L << 18) - 1L)) | (((v1 >> 18) & ((1L << 18) - 1L)) << 18);
    mTemplateCall.templateCall(endPosition, hash(l3), 3);

    final long l4 =
    (v0 & ((1L << 12) - 1L)) | ((v1 & ((1L << 12) - 1L)) << 18)
    | (((v0 >> 18) & ((1L << 6) - 1L)) << 12) | (((v1 >> 18) & ((1L << 6) - 1L)) << 30);
    mTemplateCall.templateCall(endPosition, hash(l4), 4);

    final long l5 =
    ((v0 >> 6) & ((1L << 12) - 1L)) | (((v1 >> 6) & ((1L << 12) - 1L)) << 18)
    | (((v0 >> 24) & ((1L << 6) - 1L)) << 12) | (((v1 >> 24) & ((1L << 6) - 1L)) << 30);
    mTemplateCall.templateCall(endPosition, hash(l5), 5);

    final long l6 =
    ((v0 >> 12) & ((1L << 12) - 1L)) | (((v1 >> 12) & ((1L << 12) - 1L)) << 18)
    | (((v0 >> 30) & ((1L << 6) - 1L)) << 12) | (((v1 >> 30) & ((1L << 6) - 1L)) << 30);
    mTemplateCall.templateCall(endPosition, hash(l6), 6);

    final long l7 =
    (v0 & ((1L << 12) - 1L)) | ((v1 & ((1L << 12) - 1L)) << 18)
    | (((v0 >> 24) & ((1L << 6) - 1L)) << 12) | (((v1 >> 24) & ((1L << 6) - 1L)) << 30);
    mTemplateCall.templateCall(endPosition, hash(l7), 7);

    final long l8 =
    ((v0 >> 6) & ((1L << 12) - 1L)) | (((v1 >> 6) & ((1L << 12) - 1L)) << 18)
    | (((v0 >> 30) & ((1L << 6) - 1L)) << 12) | (((v1 >> 30) & ((1L << 6) - 1L)) << 30);
    mTemplateCall.templateCall(endPosition, hash(l8), 8);

    final long l9 =
    (v0 & ((1L << 12) - 1L)) | ((v1 & ((1L << 12) - 1L)) << 18)
    | (((v0 >> 30) & ((1L << 6) - 1L)) << 12) | (((v1 >> 30) & ((1L << 6) - 1L)) << 30);
    mTemplateCall.templateCall(endPosition, hash(l9), 9);

    final long l10 =
    (v0 & ((1L << 6) - 1L)) | ((v1 & ((1L << 6) - 1L)) << 18)
    | (((v0 >> 12) & ((1L << 12) - 1L)) << 6) | (((v1 >> 12) & ((1L << 12) - 1L)) << 24);
    mTemplateCall.templateCall(endPosition, hash(l10), 10);

    final long l11 =
    ((v0 >> 6) & ((1L << 6) - 1L)) | (((v1 >> 6) & ((1L << 6) - 1L)) << 18)
    | (((v0 >> 18) & ((1L << 12) - 1L)) << 6) | (((v1 >> 18) & ((1L << 12) - 1L)) << 24);
    mTemplateCall.templateCall(endPosition, hash(l11), 11);

    final long l12 =
    ((v0 >> 12) & ((1L << 6) - 1L)) | (((v1 >> 12) & ((1L << 6) - 1L)) << 18)
    | (((v0 >> 24) & ((1L << 12) - 1L)) << 6) | (((v1 >> 24) & ((1L << 12) - 1L)) << 24);
    mTemplateCall.templateCall(endPosition, hash(l12), 12);

    final long l13 =
    (v0 & ((1L << 6) - 1L)) | ((v1 & ((1L << 6) - 1L)) << 18)
    | (((v0 >> 12) & ((1L << 6) - 1L)) << 6) | (((v1 >> 12) & ((1L << 6) - 1L)) << 24)
    | (((v0 >> 24) & ((1L << 6) - 1L)) << 12) | (((v1 >> 24) & ((1L << 6) - 1L)) << 30);
    mTemplateCall.templateCall(endPosition, hash(l13), 13);

    final long l14 =
    ((v0 >> 6) & ((1L << 6) - 1L)) | (((v1 >> 6) & ((1L << 6) - 1L)) << 18)
    | (((v0 >> 18) & ((1L << 6) - 1L)) << 6) | (((v1 >> 18) & ((1L << 6) - 1L)) << 24)
    | (((v0 >> 30) & ((1L << 6) - 1L)) << 12) | (((v1 >> 30) & ((1L << 6) - 1L)) << 30);
    mTemplateCall.templateCall(endPosition, hash(l14), 14);

    final long l15 =
    (v0 & ((1L << 6) - 1L)) | ((v1 & ((1L << 6) - 1L)) << 18)
    | (((v0 >> 12) & ((1L << 6) - 1L)) << 6) | (((v1 >> 12) & ((1L << 6) - 1L)) << 24)
    | (((v0 >> 30) & ((1L << 6) - 1L)) << 12) | (((v1 >> 30) & ((1L << 6) - 1L)) << 30);
    mTemplateCall.templateCall(endPosition, hash(l15), 15);

    final long l16 =
    (v0 & ((1L << 6) - 1L)) | ((v1 & ((1L << 6) - 1L)) << 18)
    | (((v0 >> 18) & ((1L << 12) - 1L)) << 6) | (((v1 >> 18) & ((1L << 12) - 1L)) << 24);
    mTemplateCall.templateCall(endPosition, hash(l16), 16);

    final long l17 =
    ((v0 >> 6) & ((1L << 6) - 1L)) | (((v1 >> 6) & ((1L << 6) - 1L)) << 18)
    | (((v0 >> 24) & ((1L << 12) - 1L)) << 6) | (((v1 >> 24) & ((1L << 12) - 1L)) << 24);
    mTemplateCall.templateCall(endPosition, hash(l17), 17);

    final long l18 =
    (v0 & ((1L << 6) - 1L)) | ((v1 & ((1L << 6) - 1L)) << 18)
    | (((v0 >> 18) & ((1L << 6) - 1L)) << 6) | (((v1 >> 18) & ((1L << 6) - 1L)) << 24)
    | (((v0 >> 30) & ((1L << 6) - 1L)) << 12) | (((v1 >> 30) & ((1L << 6) - 1L)) << 30);
    mTemplateCall.templateCall(endPosition, hash(l18), 18);

    final long l19 =
    (v0 & ((1L << 6) - 1L)) | ((v1 & ((1L << 6) - 1L)) << 18)
    | (((v0 >> 24) & ((1L << 12) - 1L)) << 6) | (((v1 >> 24) & ((1L << 12) - 1L)) << 24);
    mTemplateCall.templateCall(endPosition, hash(l19), 19);

    mTemplateCall.done();
  }

  @Override
  public int numberWindows() {
    return NUMBER_WINDOWS;
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("MaskL36w18s3e1 l=" + LENGTH + " w=" + NUMBER_NT + " s=3 e=1");
  }
}
