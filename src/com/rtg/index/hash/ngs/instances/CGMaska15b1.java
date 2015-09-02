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
 * <a href="doc-files/CGMaska15.jpg">This</a> illustrates the layout of the
 * masks and the calculation of the shift directions and magnitudes.
 */
public class CGMaska15b1 extends AbstractCGMask {

  private static final int READ_LENGTH = 35;

  private static final int WINDOW_LENGTH = 18;

  static final int NUMBER_BITS = 36;

  private static final int NUMBER_WINDOWS = 6;

  /**
   * A factory for the outer class.
   */
  static class CG15HashFunctionFactory implements HashFunctionFactory {
    @Override
    public NgsHashFunction create(final ReadCall readCall, final TemplateCall templateCall) {
      return new CGMaska15b1(readCall, templateCall);
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
  public static final HashFunctionFactory FACTORY = new CG15HashFunctionFactory();

  /**
   * @param readCall used in subclasses to process results of read hits.
   * @param templateCall used in subclasses to process results of template hits.
   */
  public CGMaska15b1(final ReadCall readCall, final TemplateCall templateCall) {
    super(READ_LENGTH, WINDOW_LENGTH, readCall, templateCall, LONG_BITS - (READ_LENGTH + 7));
    setHashFunction();
  }

  private static final long MASK_AR =  ((1L << 18) - 1) << 17;
  private static final long MASK_A0T =  ((1L << 5) - 1) << 37;
  private static final long MASK_A1T = ((1L << 13) - 1) << 24;
  private static final long MASK_BR =   (1L << 17) - 1;
  private static final long MASK_B0T =  ((1L << 7) - 1) << 17;
  private static final long MASK_B1T =  (1L << 10) - 1;
  private static final long MASK_CR =  ((1L << 17) - 1) << 13;
  private static final long MASK_CT =  ((1L << 17) - 1) << 20;
  private static final long MASK_DR =  ((1L << 17) - 1) << 10;
  private static final long MASK_DT =  ((1L << 17) - 1) << 17;
  private static final long MASK_E0R = ((1L <<  5) - 1) << 30;
  private static final long MASK_E1R = ((1L << 12) - 1) << 10;
  private static final long MASK_E0T = ((1L <<  5) - 1) << 37;
  private static final long MASK_E1T = ((1L << 12) - 1) << 17;
  private static final long MASK_F0R = ((1L <<  7) - 1) << 23;
  private static final long MASK_F1R =  (1L << 10) - 1;
  private static final long MASK_F0T = ((1L <<  7) - 1) << 30;
  private static final long MASK_F1T =  (1L << 10) - 1;

  @Override
  public void readAll(final int readId, long v0, long v1) {
    //System.err.println("readAll sofar=" + mSoFar);
    //System.err.println("    v0=" + com.rtg.util.Utils.toBitsSep(v0));
    //System.err.println("    v1=" + com.rtg.util.Utils.toBitsSep(v1));
    if (mSoFar < READ_LENGTH) {
      return;
    }

    // variables are named l <maskname> <mask subscript> <v0/v1>
    final long la0 = (v0 & MASK_AR) << 1;
    final long la1 = (v1 & MASK_AR) >>> 17;
    final long la = la0 | la1;
    mReadCall.readCall(readId, hash(la), 0);

    final long lb0 = (v0 & MASK_BR) << 18;
    final long lb1 =  v1 & MASK_BR;
    final long lb = lb0 | lb1;
    mReadCall.readCall(readId, hash(lb), 1);

    final long lc0 = (v0 & MASK_CR) << 5;
    final long lc1 = (v1 & MASK_CR) >>> 13;
    final long lc = lc0 | lc1;
    mReadCall.readCall(readId, hash(lc), 2);

    final long ld0 = (v0 & MASK_DR) << 8;
    final long ld1 = (v1 & MASK_DR) >>> 10;
    final long ld = ld0 | ld1;
    mReadCall.readCall(readId, hash(ld), 3);

    final long le10 = (v0 & MASK_E1R) << 8;
    final long le11 = (v1 & MASK_E1R) >>> 10;
    final long le00 =  v0 & MASK_E0R;
    final long le01 = (v1 & MASK_E0R) >>> 18;
  final long le = le00 | le10 | le01 | le11;
  mReadCall.readCall(readId, hash(le), 4);

  final long lf10 = (v0 & MASK_F1R) << 18;
  final long lf11 =  v1 & MASK_F1R;
  final long lf00 = (v0 & MASK_F0R) << 5;
  final long lf01 = (v1 & MASK_F0R) >>> 13;
  final long lf = lf00 | lf10 | lf01 | lf11;
  mReadCall.readCall(readId, hash(lf), 5);

  }

  static {
    Exam.assertEquals(MASK_AR, (MASK_A0T | MASK_A1T) >>> 7);
    Exam.assertEquals(MASK_BR, MASK_B0T >>> 7 | MASK_B1T);
  }

  @Override
  public void templateAll(final int endPosition, final long v0, final long v1) throws IOException {
    //System.err.println("templateAll sofar=" + mSoFar + " endPosition=" + endPosition);
    //System.err.println("    v0=" + com.rtg.util.Utils.toBitsSep(v0));
    //System.err.println("    v1=" + com.rtg.util.Utils.toBitsSep(v1));

    // allows for overlap of 3 and gap size of 4
    if (mSoFar < READ_LENGTH + 1) {
      return;
    }

    //The -7 allows for a maxiumum insert of 7nt the +2 is experimenntal
    final int endP = endPosition - 7 + 2;

    //System.err.println("templateAll index 0");
    // variables are named l <maskname> <mask subscript> <v0/v1> <shift x,y,z>
    //used by next 3 masks
    final long la10 = (v0 & MASK_A1T) >>> 6;
  final long la11 = (v1 & MASK_A1T) >>> 24;
  final long la = la10 | la11;

  //overlap 0
  final long la00x = (v0 & MASK_A0T) >>> 6;
  final long la01x = (v1 & MASK_A0T) >>> 24;
  final long lax = la | la00x | la01x;
  mTemplateCall.templateCall(endP, hash(lax), 0);

  //overlap -1
  final long la00y = (v0 & (MASK_A0T >>> 1)) >>> 5;
  final long la01y = (v1 & (MASK_A0T >>> 1)) >>> 23;
  final long lay = la | la00y | la01y;
  mTemplateCall.templateCall(endP/* + 1*/, hash(lay), 0);

  //overlap -2
  final long la00z = (v0 & (MASK_A0T >>> 2)) >>> 4;
  final long la01z = (v1 & (MASK_A0T >>> 2)) >>> 22;
  final long laz = la | la00z | la01z;
  mTemplateCall.templateCall(endP/* + 2*/, hash(laz), 0);

  //System.err.println("templateAll index 1");
  //used by next 3 masks
  final long lb10 = (v0 & MASK_B1T) << 18;
  final long lb11 =  v1 & MASK_B1T;
  final long lb = lb10 | lb11;

  //gap = 7
  final long lb00x = (v0 & MASK_B0T) << 11;
  final long lb01x = (v1 & MASK_B0T) >>> 7;
  final long lbx = lb | lb00x | lb01x;
  mTemplateCall.templateCall(endP/* + 1*/, hash(lbx), 1);

  //gap = 6
  final long lb00y = (v0 & (MASK_B0T >>> 1)) << 12;
  final long lb01y = (v1 & (MASK_B0T >>> 1)) >>> 6;
  final long lby = lb | lb00y | lb01y;
  mTemplateCall.templateCall(endP/* + 2*/, hash(lby), 1);

  //gap = 5
  final long lb00z = (v0 & (MASK_B0T >>> 2)) << 13;
  final long lb01z = (v1 & (MASK_B0T >>> 2)) >>> 5;
  final long lbz = lb | lb00z | lb01z;
  mTemplateCall.templateCall(endP/* + 3*/, hash(lbz), 1);

  final long lc0 = (v0 & MASK_CT) >>> 2;
  final long lc1 = (v1 & MASK_CT) >>> 20;
  final long lc = lc0 | lc1;
  mTemplateCall.templateCall(endP/* + 1*/, hash(lc), 2);

  final long ld0 = (v0 & MASK_DT) << 1;
  final long ld1 = (v1 & MASK_DT) >>> 17;
  final long ld = ld0 | ld1;
  mTemplateCall.templateCall(endP/* + 1*/, hash(ld), 3);

  //used by next 3 masks
  final long le10 = (v0 & MASK_E1T) << 1;
  final long le11 = (v1 & MASK_E1T) >>> 17;
  final long le = le10 | le11;

  //overlap = 0
  final long le00x = (v0 & MASK_E0T) >>> 7;
  final long le01x = (v1 & MASK_E0T) >>> 25;
  final long lex = le00x | le01x | le;
  mTemplateCall.templateCall(endP, hash(lex), 4);

  //overlap = -1
  final long le00y = (v0 & (MASK_E0T >>> 1)) >>> 6;
  final long le01y = (v1 & (MASK_E0T >>> 1)) >>> 24;
  final long ley = le00y | le01y | le;
  mTemplateCall.templateCall(endP/* + 1*/, hash(ley), 4);

  //overlap = -2
  final long le00z = (v0 & (MASK_E0T >>> 2)) >>> 5;
  final long le01z = (v1 & (MASK_E0T >>> 2)) >>> 23;
  final long lez = le00z | le01z | le;
  mTemplateCall.templateCall(endP/* + 2*/, hash(lez), 4);

  //used by next 3 masks
  final long lf10 = (v0 & MASK_F1T) << 18;
  final long lf11 =  v1 & MASK_F1T;
  final long lf = lf10 | lf11;

  //gap = 7
  final long lf00x = (v0 & MASK_F0T) >>> 2;
  final long lf01x = (v1 & MASK_F0T) >>> 20;
  final long lfx = lf00x | lf01x | lf;
  mTemplateCall.templateCall(endP/* + 1*/, hash(lfx), 5);

  //gap = 6
  final long lf00y = (v0 & (MASK_F0T >>> 1)) >>> 1;
  final long lf01y = (v1 & (MASK_F0T >>> 1)) >>> 19;
  final long lfy = lf00y | lf01y | lf;
  mTemplateCall.templateCall(endP/* + 2*/, hash(lfy), 5);

  //gap = 5
  final long lf00z =  v0 & (MASK_F0T >>> 2);
  final long lf01z = (v1 & (MASK_F0T >>> 2)) >>> 18;
  final long lfz = lf00z | lf01z | lf;
  mTemplateCall.templateCall(endP/* + 3*/, hash(lfz), 5);

  mTemplateCall.done();
  }

  @Override
  public int numberWindows() {
    return NUMBER_WINDOWS;
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("CG l=").append(readLength()).append(" w=").append(windowSize()).append(" a=1.5");
  }
}

