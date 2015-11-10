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
import com.rtg.util.Utils;

/**
 * r=29 w=18
 * Covers all single mismatches and most double mismatches, for overlaps of 3, 4, 5.
 * Note: There is a PDF covering these particular masks.
 */
public class CG2Maskw18a1 extends AbstractCG2Mask {

  private static final int READ_LENGTH = CgUtils.CG2_RAW_READ_LENGTH;

  private static final int WINDOW_LENGTH = 18;

  private static final int NUMBER_BITS = WINDOW_LENGTH << 1;

  private static final int NUMBER_WINDOWS = 4;

  /**
   * A factory for the outer class.
   */
  public static final HashFunctionFactory FACTORY = new HashFunctionFactory() {
    @Override
    public NgsHashFunction create(final ReadCall readCall, final TemplateCall templateCall) {
      return new CG2Maskw18a1(readCall, templateCall);
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

  private static final long MASK_A0 = 0b00000000001111111111100000000;
  private static final long MASK_A1 = 0b00000000000000000000001111111;

  private static final long MASK_B0 = 0b11111111110000000000000000000;
  private static final long MASK_B1 = 0b00000000000000000000011111111;

  private static final long MASK_C0 = 0b11111100000000000000000000000;
  private static final long MASK_C1 = 0b00000000001111111111110000000;

  private static final long MASK_D0 = 0b11111000000000000000000000000;
  private static final long MASK_D1 = 0b00000000000000001111111111111;

  /**
   * @param readCall used in subclasses to process results of read hits.
   * @param templateCall used in subclasses to process results of template hits.
   */
  public CG2Maskw18a1(final ReadCall readCall, final TemplateCall templateCall) {
    super(READ_LENGTH, WINDOW_LENGTH, readCall, templateCall, LONG_BITS - (READ_LENGTH - 4));
    setHashFunction();
    /*
    System.err.println("  MASKA=" + Utils.toBits(MASK_A, READ_LENGTH));
    System.err.println("  MASKB=" + Utils.toBits(MASK_B0 | MASK_B1, READ_LENGTH));
    System.err.println("  MASKC=" + Utils.toBits(MASK_C0 | MASK_C1, READ_LENGTH));
    */
  }

  @Override
  public void readAll(final int readId, long v0, long v1) {
    if (mSoFar < READ_LENGTH) {
      return;
    }
//    System.err.println("readAll sofar=" + mSoFar);
//    System.err.println(" v0=" + Utils.toBitsSep(v0));
//    System.err.println(" v1=" + Utils.toBitsSep(v1));
    final long laa0 = (v0 & MASK_A0) << 17;
    final long laa1 = (v0 & MASK_A1) << 18;
    final long lba0 = (v1 & MASK_A0) >> 1;
    final long lba1 =  v1 & MASK_A1;
    final long packedA = laa0 | lba0 | laa1 | lba1;
    mReadCall.readCall(readId, hash(packedA), 0);

    final long lab0 = (v0 & MASK_B0) << 7;
    final long lab1 = (v0 & MASK_B1) << 18;
    final long lbb0 = (v1 & MASK_B0) >> 11;
    final long lbb1 = v1 & MASK_B1;
    final long packedB = lab0 | lbb0 | lab1 | lbb1;
    mReadCall.readCall(readId, hash(packedB), 1);

    final long lac0 = (v0 & MASK_C0) << 7;
    final long lac1 = (v0 & MASK_C1) << 11;
    final long lbc0 = (v1 & MASK_C0) >> 11;
    final long lbc1 = (v1 & MASK_C1) >> 7;
    final long packedC = lac0 | lbc0 | lac1 | lbc1;
    mReadCall.readCall(readId, hash(packedC), 2);

    final long lad0 = (v0 & MASK_D0) << 7;
    final long lad1 = (v0 & MASK_D1) << 18;
    final long lbd0 = (v1 & MASK_D0) >> 11;
    final long lbd1 = v1 & MASK_D1;
    final long packedD = lad0 | lbd0 | lad1 | lbd1;
    mReadCall.readCall(readId, hash(packedD), 3);

    //dumpBlocks("RMASKD", lad0, lad1, lbd0, lbd1);
  }

  private static final long MASKT_A0 = MASK_A0;
  private static final long MASKT_A1 = MASK_A1;

  private static final long MASKT_B0_3 = MASK_B0 >> 3;
  private static final long MASKT_B0_4 = MASKT_B0_3 >> 1;
  private static final long MASKT_B0_5 = MASKT_B0_4 >> 1;
  private static final long MASKT_B1 = MASK_B1;

  private static final long MASKT_C0_3 = MASK_C0 >> 3;
  private static final long MASKT_C0_4 = MASKT_C0_3 >> 1;
  private static final long MASKT_C0_5 = MASKT_C0_4 >> 1;
  private static final long MASKT_C1 = MASK_C1;

  private static final long MASKT_D0_3 = MASK_D0 >> 3;
  private static final long MASKT_D0_4 = MASKT_D0_3 >> 1;
  private static final long MASKT_D0_5 = MASKT_D0_4 >> 1;
  private static final long MASKT_D1 = MASK_D1;


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

    // index 0, doesn't care about overlap variation
    final long laa0 = (v0 & MASKT_A0) << 17;
    final long laa1 = (v0 & MASKT_A1) << 18;
    final long lba0 = (v1 & MASKT_A0) >> 1;
    final long lba1 =  v1 & MASKT_A1;
    final long packedA = laa0 | lba0 | laa1 | lba1;
    mTemplateCall.templateCall(endP, hash(packedA), 0);


    // index 1, most common overlaps
    // Overlap of 3
    final long lab03 = (v0 & MASKT_B0_3) << 10;
    final long lbb03 = (v1 & MASKT_B0_3) >> 8;
    final long lab1 = (v0 & MASKT_B1) << 18;
    final long lbb1 = v1 & MASKT_B1;
    final long packedB3 = lab03 | lbb03 | lab1 | lbb1;
    mTemplateCall.templateCall(endP, hash(packedB3), 1);

    // Overlap of 4
    final long lab04 = (v0 & MASKT_B0_4) << 11;
    final long lbb04 = (v1 & MASKT_B0_4) >> 7;
    final long packedB4 = lab04 | lbb04 | lab1 | lbb1;
    mTemplateCall.templateCall(endP, hash(packedB4), 1);

    // Overlap of 5
    final long lab05 = (v0 & MASKT_B0_5) << 12;
    final long lbb05 = (v1 & MASKT_B0_5) >> 6;
    final long packedB5 = lab05 | lbb05 | lab1 | lbb1;
    mTemplateCall.templateCall(endP, hash(packedB5), 1);

    // index 2, most common overlaps
    // Overlap of 3
    final long lac03 = (v0 & MASKT_C0_3) << 10;
    final long lbc03 = (v1 & MASKT_C0_3) >> 8;
    final long lac1 = (v0 & MASKT_C1) << 11;
    final long lbc1 = (v1 & MASKT_C1) >> 7;
    final long packedC3 = lac03 | lbc03 | lac1 | lbc1;
    mTemplateCall.templateCall(endP, hash(packedC3), 2);

    // Overlap of 4
    final long lac04 = (v0 & MASKT_C0_4) << 11;
    final long lbc04 = (v1 & MASKT_C0_4) >> 7;
    final long packedC4 = lac04 | lbc04 | lac1 | lbc1;
    mTemplateCall.templateCall(endP, hash(packedC4), 2);

    // Overlap of 5
    final long lac05 = (v0 & MASKT_C0_5) << 12;
    final long lbc05 = (v1 & MASKT_C0_5) >> 6;
    final long packedC5 = lac05 | lbc05 | lac1 | lbc1;
    mTemplateCall.templateCall(endP, hash(packedC5), 2);

    // index 3, most common overlaps
    // Overlap of 3
    final long lad03 = (v0 & MASKT_D0_3) << 10;
    final long lbd03 = (v1 & MASKT_D0_3) >> 8;
    final long lad1 = (v0 & MASKT_D1) << 18;
    final long lbd1 = v1 & MASKT_D1;
    final long packedD3 = lad03 | lbd03 | lad1 | lbd1;
    mTemplateCall.templateCall(endP, hash(packedD3), 3);

    // Overlap of 4
    final long lad04 = (v0 & MASKT_D0_4) << 11;
    final long lbd04 = (v1 & MASKT_D0_4) >> 7;
    final long packedD4 = lad04 | lbd04 | lad1 | lbd1;
    mTemplateCall.templateCall(endP, hash(packedD4), 3);

    // Overlap of 5
    final long lad05 = (v0 & MASKT_D0_5) << 12;
    final long lbd05 = (v1 & MASKT_D0_5) >> 6;
    final long packedD5 = lad05 | lbd05 | lad1 | lbd1;
    mTemplateCall.templateCall(endP, hash(packedD5), 3);

    //dumpBlocks("TMASKA", laa0, laa1, lba0, lba1);
    //dumpBlocks("TMASKD", lad03, lad1, lbd03, lbd1);

    mTemplateCall.done();
  }

  protected static void dumpBlocks(String label, long... blocks) {
    for (int i = 0; i < blocks.length; i++) {
      System.err.println(label + i + "=" + Utils.toBits(blocks[i], NUMBER_BITS));
    }
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

