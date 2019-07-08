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
import com.rtg.util.integrity.Exam;

/**
 * r=35 s=1
 * Allows for 0, -1, -2 on left hand gap and 5, 6, 7 on right hand gap.
 * Misses % of possibilities at varying substitutions:
 * <table><caption>fraction missed by substitutions</caption>
 * <tr> <td>1</td> <td>0.0%</td> </tr>
 * <tr> <td>2</td> <td>51.4%</td> </tr>
 * <tr> <td>3</td> <td>77.1%</td> </tr>
 * <tr> <td>3</td> <td>89.6%</td> </tr>
 * </table>
 *
 */
public class CGMaska1b1 extends AbstractCGMask {

  private static final int READ_LENGTH = 35;

  private static final int WINDOW_LENGTH = 18;

  static final int NUMBER_BITS = 36;

  private static final int NUMBER_WINDOWS = 2;

  /**
   * A factory for the outer class.
   */
  static class CG1HashFunctionFactory implements CGHashFunctionFactory {
    @Override
    public NgsHashFunction create(final ReadCall readCall, final TemplateCall templateCall) {
      return new CGMaska1b1(readCall, templateCall);
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
  }

  /**
   * A factory for the outer class.
   */
  public static final HashFunctionFactory FACTORY = new CG1HashFunctionFactory();

  /**
   * @param readCall used in subclasses to process results of read hits.
   * @param templateCall used in subclasses to process results of template hits.
   */
  public CGMaska1b1(final ReadCall readCall, final TemplateCall templateCall) {
    super(WINDOW_LENGTH, readCall, templateCall);
    setHashFunction();
  }

  private static final long MASK_A = ((1L << 18) - 1) << (7 + 10);
  private static final long MASK_B = (1L << 17) - 1;

  @Override
  public void readAll(final int readId, long v0, long v1) {
    //System.err.println("readAll sofar=" + mSoFar);
    //System.err.println("    v0=" + com.rtg.util.Utils.toBitsSep(v0));
    //System.err.println("    v1=" + com.rtg.util.Utils.toBitsSep(v1));
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

  /** Masks for left 18 nt of read (includes the 5 nt region and 0, -1, -2 gap). */
  private static final long MASKT_A0 = ((1L << 5) - 1) << (20 + 7 + 10);
  private static final long MASKT_A1 = ((1L << 13) - 1) << (7 + 7 + 10);
  static {
    Exam.assertEquals((MASKT_A0 | MASKT_A1) >>> 7, MASK_A);
    Exam.assertEquals(MASKT_A0 & MASKT_A1, 0);
    Exam.assertTrue((MASKT_A0 >>> 1 & MASKT_A1) != 0);
  }

  /** Masks for right 17 nt of read (includes the 10 nt region and 5, 6, 7 gap). */
  private static final long MASKT_B0 = (1L << 10) - 1;
  private static final long MASKT_B1 = ((1L << 7) - 1) << 17;

  @Override
  public void templateAll(final int endPosition, final long v0, final long v1) throws IOException {
    //System.err.println("templateAll sofar=" + mSoFar + " endPosition=" + endPosition);
    //System.err.println("    v0=" + com.rtg.util.Utils.toBitsSep(v0));
    //System.err.println("    v1=" + com.rtg.util.Utils.toBitsSep(v1));
    if (mSoFar < READ_LENGTH + 1) {
      return;
    }

    //The -7 allows for a maxiumum insert of 7nt the +2 is experimenntal
    final int endP = endPosition - 7 + 2;

    //System.err.println("templateAll index 0");
    //lb used by next 3 masks
    final long la10 = (v0 & MASKT_A1) >>> 6;
  final long la11 = (v1 & MASKT_A1) >>> 24;
  final long la = la10 | la11;

  //overlap 0
  final long la00x = (v0 & MASKT_A0) >>> 6;
  final long la01x = (v1 & MASKT_A0) >>> 24;
    final long lax = la | la00x | la01x;
    mTemplateCall.templateCall(endP, hash(lax), 0);

    //overlap -1
    final long la00y = (v0 & (MASKT_A0 >>> 1)) >>> 5;
  final long la01y = (v1 & (MASKT_A0 >>> 1)) >>> 23;
  final long lay = la | la00y | la01y;
  mTemplateCall.templateCall(endP/* + 1*/, hash(lay), 0);

  //overlap -2
  final long la00z = (v0 & (MASKT_A0 >>> 2)) >>> 4;
  final long la01z = (v1 & (MASKT_A0 >>> 2)) >>> 22;
  final long laz = la | la00z | la01z;
  mTemplateCall.templateCall(endP/* + 2*/, hash(laz), 0);

  //System.err.println("templateAll index 1");
  //ld used by next 3 masks
  final long lb10 = (v0 & MASKT_B0) << 18;
  final long lb11 =  v1 & MASKT_B0;
  final long lb = lb10 | lb11;

  //gap = 7
  final long lb00x = (v0 & MASKT_B1) << 11;
  final long lb01x = (v1 & MASKT_B1) >>> 7;
  final long lbx = lb | lb00x | lb01x;
  mTemplateCall.templateCall(endP/* + 1*/, hash(lbx), 1);

  //gap = 6
  final long lb00y = (v0 & (MASKT_B1 >>> 1)) << 12;
  final long lb01y = (v1 & (MASKT_B1 >>> 1)) >>> 6;
  final long lby = lb | lb00y | lb01y;
  mTemplateCall.templateCall(endP/* + 2*/, hash(lby), 1);

  //gap = 5
  final long lb00z = (v0 & (MASKT_B1 >>> 2)) << 13;
  final long lb01z = (v1 & (MASKT_B1 >>> 2)) >>> 5;
  final long lbz = lb | lb00z | lb01z;
  mTemplateCall.templateCall(endP/* + 3*/, hash(lbz), 1);

  mTemplateCall.done();
  }

  @Override
  public int numberWindows() {
    return NUMBER_WINDOWS;
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("CG l=").append(readLength()).append(" w=").append(windowSize()).append(" e=1");
  }
}

