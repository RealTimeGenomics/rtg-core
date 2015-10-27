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
import com.rtg.reader.CgUtils;

/**
 * r=29 s=0
 * Covers all single mismatches, for overlaps of 3, 4, 5.
 */
public class CG2Maska15 extends AbstractCGMask {

  private static final int READ_LENGTH = CgUtils.CG2_RAW_READ_LENGTH;

  private static final int WINDOW_LENGTH = 15;

  private static final int NUMBER_BITS = 30;

  private static final int NUMBER_WINDOWS = 5;

  /**
   * A factory for the outer class.
   */
  public static final HashFunctionFactory FACTORY = new HashFunctionFactory() {
    @Override
    public NgsHashFunction create(final ReadCall readCall, final TemplateCall templateCall) {
      return new CG2Maska15(readCall, templateCall);
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

  private static final long MASK_A0 = 0b11111111110000000000000000000;
  private static final long MASK_A1 = 0b00000000000000000000000011111;
  private static final long MASK_B  = 0b00000000001111111111111100000;
  private static final long MASK_C  = 0b00000000000000011111111111111;
  private static final long MASK_D0 = 0b00000111110000000000000000000;
  private static final long MASK_D1 = 0b00000000000000011111111100000;
  private static final long MASK_E0 = 0b11111000000000000000000000000;
  private static final long MASK_E1 = 0b00000000000000011111000000000;
  private static final long MASK_E2 = 0b00000000000000000000000001111;

  /**
   * @param readCall used in subclasses to process results of read hits.
   * @param templateCall used in subclasses to process results of template hits.
   */
  public CG2Maska15(final ReadCall readCall, final TemplateCall templateCall) {
    super(READ_LENGTH, WINDOW_LENGTH, readCall, templateCall, LONG_BITS - (READ_LENGTH - 4));
    setHashFunction();

//    System.err.println(" MASKA=" + Utils.toBits(MASK_A0 | MASK_A1, READ_LENGTH));
//    System.err.println(" MASKB=" + Utils.toBits(MASK_B, READ_LENGTH));
//    System.err.println(" MASKC=" + Utils.toBits(MASK_C, READ_LENGTH));
//    System.err.println(" MASKD=" + Utils.toBits(MASK_D0 | MASK_D1, READ_LENGTH));
//    System.err.println(" MASKE=" + Utils.toBits(MASK_E0 | MASK_E1 | MASK_E2, READ_LENGTH));
//
//    System.err.println("MASKTA=" + Utils.toBits(MASKT_A0_1 | MASKT_A1, READ_LENGTH));
//    System.err.println("MASKTB=" + Utils.toBits(MASKT_B, READ_LENGTH));
//    System.err.println("MASKTC=" + Utils.toBits(MASKT_C, READ_LENGTH));
//    System.err.println("MASKTD=" + Utils.toBits(MASKT_D0_1 | MASKT_D1, READ_LENGTH));
//    System.err.println("MASKTE=" + Utils.toBits(MASKT_E0_1 | MASKT_E1 | MASKT_E2, READ_LENGTH));
  }

  @Override
  public void readAll(final int readId, long v0, long v1) {
    if (mSoFar < READ_LENGTH) {
      return;
    }
//    System.err.println("readAll sofar=" + mSoFar);
//    System.err.println(" v0=" + Utils.toBitsSep(v0));
//    System.err.println(" v1=" + Utils.toBitsSep(v1));
    final long la0 = (v0 & MASK_A0) << 1;
    final long la1 = (v0 & MASK_A1) << 15;
    final long lb0 = (v1 & MASK_A0) >> 14;
    final long lb1 = v1 & MASK_A1;
    final long l2x0 = la0 | lb0 | la1 | lb1;
    mReadCall.readCall(readId, hash(l2x0), 0);

    final long lc = (v0 & MASK_B) << (14 - 5);
    final long ld =  (v1 & MASK_B) >> 5;
    final long l2x1 = lc | ld;
    mReadCall.readCall(readId, hash(l2x1), 1);

    final long le = (v0 & MASK_C) << 14;
    final long lf = v1 & MASK_C;
    final long l2x2 = le | lf;
    mReadCall.readCall(readId, hash(l2x2), 2);

    final long lg0 = (v0 & MASK_D0) << 4;
    final long lg1 = (v0 & MASK_D1) << 9;
    final long lh0 = (v1 & MASK_D0) >> (19 - 9);
    final long lh1 = (v1 & MASK_D1) >> 5;
    final long l2x3 = lg0 | lh0 | lg1 | lh1;
    mReadCall.readCall(readId, hash(l2x3), 3);

    final long li0 = (v0 & MASK_E0) >> 1;
    final long li1 = (v0 & MASK_E1) << 9;
    final long li2 = (v0 & MASK_E2) << 14;
    final long lj0 = (v1 & MASK_E0) >> 15;
    final long lj1 = (v1 & MASK_E1) >> 5;
    final long lj2 = v1 & MASK_E2;
    final long l2x4 = li0 | lj0 | li1 | lj1 | li2 | lj2;
    mReadCall.readCall(readId, hash(l2x4), 4);

//    System.err.println(" MASKEa=" + Utils.toBits(li0, READ_LENGTH));
//    System.err.println(" MASKEb=" + Utils.toBits(li1, READ_LENGTH));
//    System.err.println(" MASKEc=" + Utils.toBits(li2, READ_LENGTH));
//    System.err.println(" MASKEd=" + Utils.toBits(lj0, READ_LENGTH));
//    System.err.println(" MASKEe=" + Utils.toBits(lj1, READ_LENGTH));
//    System.err.println(" MASKEf=" + Utils.toBits(lj2, READ_LENGTH));
  }

  private static final long MASKT_A0_1 = MASK_A0 >> 3;
  private static final long MASKT_A0_2 = MASKT_A0_1 >> 1;
  private static final long MASKT_A0_3 = MASKT_A0_2 >> 1;
  private static final long MASKT_A1 = MASK_A1;
  private static final long MASKT_B = MASK_B;
  private static final long MASKT_C = MASK_C;
  private static final long MASKT_D0_1 = MASK_D0 >> 3;
  private static final long MASKT_D0_2 = MASKT_D0_1 >> 1;
  private static final long MASKT_D0_3 = MASKT_D0_2 >> 1;
  private static final long MASKT_D1 = MASK_D1;
  private static final long MASKT_E0_1 = MASK_E0 >> 3;
  private static final long MASKT_E0_2 = MASKT_E0_1 >> 1;
  private static final long MASKT_E0_3 = MASKT_E0_2 >> 1;
  private static final long MASKT_E1 = MASK_E1;
  private static final long MASKT_E2 = MASK_E2;

  @Override
  public void templateAll(final int endPosition, long v0, long v1) throws IOException {
//    System.err.println("templateAll sofar=" + mSoFar + " endPosition=" + endPosition);
//    System.err.println(" v0=" + Utils.toBitsSep(v0));
//    System.err.println(" v1=" + Utils.toBitsSep(v1));
    if (mSoFar < READ_LENGTH + 1) {
      return;
    }

    //The +3 adjusts for typical overlap of 3 bases (not critical)
    final int endP = endPosition + 3;

    // index 0, most common overlaps
    // Overlap of 3
    final long le0 = (v0 & MASKT_A0_1) << 4;
    final long lf0 = (v1 & MASKT_A0_1) >> 11;
    final long le1 = (v0 & MASKT_A1) << 15;
    final long lf1 = v1 & MASKT_A1;
    final long l3x0 = le0 | lf0 | le1 | lf1;
    mTemplateCall.templateCall(endP, hash(l3x0), 0);

    // Overlap of 4
    final long la0 = (v0 & MASKT_A0_2) << 5;
    final long lb0 = (v1 & MASKT_A0_2) >> 10;
    final long l4x0 = la0 | lb0 | le1 | lf1;
    mTemplateCall.templateCall(endP, hash(l4x0), 0);

    // Overlap of 5
    final long lg0 = (v0 & MASKT_A0_3) << 6;
    final long lh0 = (v1 & MASKT_A0_3) >> 9;
    final long l5x0 = lg0 | lh0 | le1 | lf1;
    mTemplateCall.templateCall(endP, hash(l5x0), 0);


    // index 1 (doesn't care about overlap variation)
    final long lc = (v0 & MASKT_B) << (14 - 5);
    final long ld =  (v1 & MASKT_B) >> 5;
    final long l0x1 = lc | ld;
    mTemplateCall.templateCall(endP, hash(l0x1), 1);


    // index 2 (doesn't care about overlap variation)
    final long li = (v0 & MASKT_C) << 14;
    final long lj = v1 & MASKT_C;
    final long l0x2 = li | lj;
    mTemplateCall.templateCall(endP, hash(l0x2), 2);


    // index 3 (most common overlaps)
    // Overlap of 3
    final long lk0 = (v0 & MASKT_D0_1) << 7;
    final long ll0 = (v1 & MASKT_D0_1) >> 7;
    final long lk1 = (v0 & MASKT_D1) << 9;
    final long ll1 = (v1 & MASKT_D1) >> 5;
    final long l3x3 = lk0 | ll0 | lk1 | ll1;
    mTemplateCall.templateCall(endP, hash(l3x3), 3);

    // Overlap of 4
    final long lm0 = (v0 & MASKT_D0_2) << 8;
    final long ln0 = (v1 & MASKT_D0_2) >> 6;
    final long l4x3 = lm0 | ln0 | lk1 | ll1;
    mTemplateCall.templateCall(endP, hash(l4x3), 3);

    // Overlap of 5
    final long lo0 = (v0 & MASKT_D0_3) << 9;
    final long lp0 = (v1 & MASKT_D0_3) >> 5;
    final long l5x3 = lo0 | lp0 | lk1 | ll1;
    mTemplateCall.templateCall(endP, hash(l5x3), 3);

    // index 4 (most common overlaps)
    // Overlap of 3
    final long lq0 = (v0 & MASKT_E0_1) << 2;
    final long lr0 = (v1 & MASKT_E0_1) >> 12;
    final long lq1 = (v0 & MASKT_E1) << 9;
    final long lr1 = (v1 & MASKT_E1) >> 5;
    final long lq2 = (v0 & MASKT_E2) << 14;
    final long lr2 = v1 & MASKT_E2;
    final long l12 = lq1 | lq2 | lr1 | lr2;
    final long l3x4 = lq0 | lr0 | l12;
    mTemplateCall.templateCall(endP, hash(l3x4), 4);

    // Overlap of 4
    final long ls0 = (v0 & MASKT_E0_2) << 3;
    final long lt0 = (v1 & MASKT_E0_2) >> 11;
    final long l4x4 = ls0 | lt0 | l12;
    mTemplateCall.templateCall(endP, hash(l4x4), 4);

    // Overlap of 5
    final long lu0 = (v0 & MASKT_E0_3) << 4;
    final long lv0 = (v1 & MASKT_E0_3) >> 10;
    final long l4x5 = lu0 | lv0 | l12;
    mTemplateCall.templateCall(endP, hash(l4x5), 4);

//    System.err.println("MASKTDa=" + Utils.toBits(lu0, READ_LENGTH));
//    System.err.println("MASKTDb=" + Utils.toBits(lq1, READ_LENGTH));
//    System.err.println("MASKTDc=" + Utils.toBits(lq2, READ_LENGTH));
//    System.err.println("MASKTDd=" + Utils.toBits(lv0, READ_LENGTH));
//    System.err.println("MASKTDe=" + Utils.toBits(lr1, READ_LENGTH));
//    System.err.println("MASKTDf=" + Utils.toBits(lr2, READ_LENGTH));

    mTemplateCall.done();
  }

  @Override
  public int numberWindows() {
    return NUMBER_WINDOWS;
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("CG2 l=").append(readLength()).append(" w=").append(windowSize()).append(" e=0");
  }
}

