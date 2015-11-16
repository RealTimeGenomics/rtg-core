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
 * Covers all single mismatches and most double mismatches, for overlaps of 3, 4, 5.
 * Note: There is a PDF covering these particular masks.
 */
public class CG2Maska15 extends AbstractCG2Mask {

  private static final int READ_LENGTH = CgUtils.CG2_RAW_READ_LENGTH;

  private static final int WINDOW_LENGTH = 15;

  private static final int NUMBER_BITS = 30;

  private static final int NUMBER_WINDOWS = 7;

  /**
   * A factory for the outer class.
   */
  public static final HashFunctionFactory FACTORY = new CG2HashFunctionFactory() {
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
  private static final long MASK_F0 = 0b10101010100000000000000000000;
  private static final long MASK_F1 = 0b00000000000010101010101010101;
  private static final long MASK_G0 = 0b01010101010000000000000000000;
  private static final long MASK_G1 = 0b00000000000101010101010101010;

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
    final long packedA = la0 | lb0 | la1 | lb1;
    mReadCall.readCall(readId, hash(packedA), 0);

    final long lc = (v0 & MASK_B) << 11;
    final long ld =  (v1 & MASK_B) >> 3;
    final long packedB = lc | ld;
    mReadCall.readCall(readId, hash(packedB), 1);

    final long le = (v0 & MASK_C) << 16;
    final long lf = (v1 & MASK_C) << 2;
    final long packedC = le | lf;
    mReadCall.readCall(readId, hash(packedC), 2);

    final long lg0 = (v0 & MASK_D0) << 6;
    final long lg1 = (v0 & MASK_D1) << 11;
    final long lh0 = (v1 & MASK_D0) >> 8;
    final long lh1 = (v1 & MASK_D1) >> 3;
    final long packedD = lg0 | lh0 | lg1 | lh1;
    mReadCall.readCall(readId, hash(packedD), 3);

    final long li0 = (v0 & MASK_E0) << 1;
    final long li1 = (v0 & MASK_E1) << 11;
    final long li2 = (v0 & MASK_E2) << 16;
    final long lj0 = (v1 & MASK_E0) >> 13;
    final long lj1 = (v1 & MASK_E1) >> 3;
    final long lj2 = (v1 & MASK_E2) << 2;
    final long packedE = li0 | lj0 | li1 | lj1 | li2 | lj2;
    mReadCall.readCall(readId, hash(packedE), 4);

    final long lk0 = (v0 & MASK_F0) << 1;
    final long ll0 = v1 & MASK_F0;
    final long lk1 = (v0 & MASK_F1) << 3;
    final long ll1 = (v1 & MASK_F1) << 2;
    final long packedF = lk0 | ll0 | lk1 | ll1;
    mReadCall.readCall(readId, hash(packedF), 5);

    final long lm0 = (v0 & MASK_G0) << 2;
    final long ln0 = (v1 & MASK_G0) << 1;
    final long lm1 = (v0 & MASK_G1) << 2;
    final long ln1 = (v1 & MASK_G1) << 1;
    final long packedG = lm0 | ln0 | lm1 | ln1;
    mReadCall.readCall(readId, hash(packedG), 6);

//    System.err.println(" MASKGa=" + Utils.toBits(lm0, READ_LENGTH));
//    System.err.println(" MASKGb=" + Utils.toBits(ln0, READ_LENGTH));
//    System.err.println(" MASKGc=" + Utils.toBits(lm1, READ_LENGTH));
//    System.err.println(" MASKGd=" + Utils.toBits(ln1, READ_LENGTH));
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
  private static final long MASKT_F0_1 = MASK_F0 >> 3;
  private static final long MASKT_F0_2 = MASKT_F0_1 >> 1;
  private static final long MASKT_F0_3 = MASKT_F0_2 >> 1;
  private static final long MASKT_F1 = MASK_F1;
  private static final long MASKT_G0_1 = MASK_G0 >> 3;
  private static final long MASKT_G0_2 = MASKT_G0_1 >> 1;
  private static final long MASKT_G0_3 = MASKT_G0_2 >> 1;
  private static final long MASKT_G1 = MASK_G1;

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
    final long packedA3 = le0 | lf0 | le1 | lf1;
    mTemplateCall.templateCall(endP, hash(packedA3), 0);

    // Overlap of 4
    final long la0 = (v0 & MASKT_A0_2) << 5;
    final long lb0 = (v1 & MASKT_A0_2) >> 10;
    final long packedA4 = la0 | lb0 | le1 | lf1;
    mTemplateCall.templateCall(endP, hash(packedA4), 0);

    // Overlap of 5
    final long lg0 = (v0 & MASKT_A0_3) << 6;
    final long lh0 = (v1 & MASKT_A0_3) >> 9;
    final long packedA5 = lg0 | lh0 | le1 | lf1;
    mTemplateCall.templateCall(endP, hash(packedA5), 0);


    // index 1 (doesn't care about overlap variation)
    final long lc = (v0 & MASKT_B) << 11;
    final long ld =  (v1 & MASKT_B) >> 3;
    final long packedB = lc | ld;
    mTemplateCall.templateCall(endP, hash(packedB), 1);


    // index 2 (doesn't care about overlap variation)
    final long li = (v0 & MASKT_C) << 16;
    final long lj = (v1 & MASKT_C) << 2;
    final long packedC = li | lj;
    mTemplateCall.templateCall(endP, hash(packedC), 2);


    // index 3 (most common overlaps)
    // Overlap of 3
    final long lk0 = (v0 & MASKT_D0_1) << 9;
    final long ll0 = (v1 & MASKT_D0_1) >> 5;
    final long lk1 = (v0 & MASKT_D1) << 11;
    final long ll1 = (v1 & MASKT_D1) >> 3;
    final long packedD3 = lk0 | ll0 | lk1 | ll1;
    mTemplateCall.templateCall(endP, hash(packedD3), 3);

    // Overlap of 4
    final long lm0 = (v0 & MASKT_D0_2) << 10;
    final long ln0 = (v1 & MASKT_D0_2) >> 4;
    final long packedD4 = lm0 | ln0 | lk1 | ll1;
    mTemplateCall.templateCall(endP, hash(packedD4), 3);

    // Overlap of 5
    final long lo0 = (v0 & MASKT_D0_3) << 11;
    final long lp0 = (v1 & MASKT_D0_3) >> 3;
    final long packedD5 = lo0 | lp0 | lk1 | ll1;
    mTemplateCall.templateCall(endP, hash(packedD5), 3);

    // index 4 (most common overlaps)
    // Overlap of 3
    final long lq0 = (v0 & MASKT_E0_1) << 4;
    final long lr0 = (v1 & MASKT_E0_1) >> 10;
    final long lq1 = (v0 & MASKT_E1) << 11;
    final long lr1 = (v1 & MASKT_E1) >> 3;
    final long lq2 = (v0 & MASKT_E2) << 16;
    final long lr2 = (v1 & MASKT_E2) << 2;
    final long l12 = lq1 | lq2 | lr1 | lr2;
    final long packedE3 = lq0 | lr0 | l12;
    mTemplateCall.templateCall(endP, hash(packedE3), 4);

    // Overlap of 4
    final long ls0 = (v0 & MASKT_E0_2) << 5;
    final long lt0 = (v1 & MASKT_E0_2) >> 9;
    final long packedE4 = ls0 | lt0 | l12;
    mTemplateCall.templateCall(endP, hash(packedE4), 4);

    // Overlap of 5
    final long lu0 = (v0 & MASKT_E0_3) << 6;
    final long lv0 = (v1 & MASKT_E0_3) >> 8;
    final long packedE5 = lu0 | lv0 | l12;
    mTemplateCall.templateCall(endP, hash(packedE5), 4);

    // index 5 (most common overlaps)
    // Overlap of 3
    final long lw0 = (v0 & MASKT_F0_1) << 4;
    final long lx0 = (v1 & MASKT_F0_1) << 3;
    final long lw1 = (v0 & MASKT_F1) << 3;
    final long lx1 = (v1 & MASKT_F1) << 2;
    final long packedF3 = lw0 | lx0 | lw1 | lx1;
    mTemplateCall.templateCall(endP, hash(packedF3), 5);

    // Overlap of 4
    final long ly0 = (v0 & MASKT_F0_2) << 5;
    final long lz0 = (v1 & MASKT_F0_2) << 4;
    final long packedF4 = ly0 | lz0 | lw1 | lx1;
    mTemplateCall.templateCall(endP, hash(packedF4), 5);

    // Overlap of 5
    final long laa0 = (v0 & MASKT_F0_3) << 6;
    final long lab0 = (v1 & MASKT_F0_3) << 7;
    final long packedF5 = laa0 | lab0 | lw1 | lx1;
    mTemplateCall.templateCall(endP, hash(packedF5), 5);

    // index 6 (most common overlaps)
    // Overlap of 3
    final long lac0 = (v0 & MASKT_G0_1) << 5;
    final long lad0 = (v1 & MASKT_G0_1) << 4;
    final long lac1 = (v0 & MASKT_G1) << 2;
    final long lad1 = (v1 & MASKT_G1) << 1;
    final long packedG3 = lac0 | lad0 | lac1 | lad1;
    mTemplateCall.templateCall(endP, hash(packedG3), 6);

    // Overlap of 4
    final long lae0 = (v0 & MASKT_G0_2) << 6;
    final long laf0 = (v1 & MASKT_G0_2) << 5;
    final long packedG4 = lae0 | laf0 | lac1 | lad1;
    mTemplateCall.templateCall(endP, hash(packedG4), 6);

    // Overlap of 5
    final long lag0 = (v0 & MASKT_G0_3) << 7;
    final long lah0 = (v1 & MASKT_G0_3) << 6;
    final long packedG5 = lag0 | lah0 | lac1 | lad1;
    mTemplateCall.templateCall(endP, hash(packedG5), 6);

//    System.err.println("MASKTGa=" + Utils.toBits(lac0, READ_LENGTH));
//    System.err.println("MASKTGb=" + Utils.toBits(lad0, READ_LENGTH));
//    System.err.println("MASKTGc=" + Utils.toBits(lac1, READ_LENGTH));
//    System.err.println("MASKTGd=" + Utils.toBits(lad1, READ_LENGTH));

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

