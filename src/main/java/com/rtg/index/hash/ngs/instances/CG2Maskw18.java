/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
package com.rtg.index.hash.ngs.instances;

import java.io.IOException;

import com.rtg.index.hash.ngs.HashFunctionFactory;
import com.rtg.index.hash.ngs.NgsHashFunction;
import com.rtg.index.hash.ngs.ReadCall;
import com.rtg.index.hash.ngs.TemplateCall;
import com.rtg.reader.CgUtils;

/**
 * r=29 w=18
 * Covers all single mismatches and most double mismatches, for overlaps of 3, 4, 5.
 * Note: There is a PDF covering these particular masks.
 */
public class CG2Maskw18 extends AbstractCG2Mask {

  private static final int READ_LENGTH = CgUtils.CG2_RAW_READ_LENGTH;

  private static final int WINDOW_LENGTH = 18;

  private static final int NUMBER_BITS = WINDOW_LENGTH << 1;

  private static final int NUMBER_WINDOWS = 3;

  /**
   * A factory for the outer class.
   */
  public static final HashFunctionFactory FACTORY = new CG2HashFunctionFactory() {
    @Override
    public NgsHashFunction create(final ReadCall readCall, final TemplateCall templateCall) {
      return new CG2Maskw18(readCall, templateCall);
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

  private static final long MASK_A  = 0b00000000000111111111111111111;

  private static final long MASK_B0 = 0b11111111110000000000000000000;
  private static final long MASK_B1 = 0b00000000000000000000011111111;

  private static final long MASK_C0 = 0b11111100000000000000000000000;
  private static final long MASK_C1 = 0b00000000001111111111110000000;

  /**
   * @param readCall used in subclasses to process results of read hits.
   * @param templateCall used in subclasses to process results of template hits.
   */
  public CG2Maskw18(final ReadCall readCall, final TemplateCall templateCall) {
    super(WINDOW_LENGTH, readCall, templateCall);
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
    final long laa = (v0 & MASK_A) << 18;
    final long lba =  v1 & MASK_A;
    final long packedA = laa | lba;
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


//    System.err.println(" MASKAa=" + Utils.toBits(laa, NUMBER_BITS));
//    System.err.println(" MASKAb=" + Utils.toBits(lba, NUMBER_BITS));

//    System.err.println(" MASKCa=" + Utils.toBits(lac0, NUMBER_BITS));
//    System.err.println(" MASKCb=" + Utils.toBits(lac1, NUMBER_BITS));
//    System.err.println(" MASKCc=" + Utils.toBits(lbc0, NUMBER_BITS));
//    System.err.println(" MASKCd=" + Utils.toBits(lbc1, NUMBER_BITS));

  }

  private static final long MASKT_A = MASK_A;

  private static final long MASKT_B0_3 = MASK_B0 >> 3;
  private static final long MASKT_B0_4 = MASKT_B0_3 >> 1;
  private static final long MASKT_B0_5 = MASKT_B0_4 >> 1;
  private static final long MASKT_B1 = MASK_B1;

  private static final long MASKT_C0_3 = MASK_C0 >> 3;
  private static final long MASKT_C0_4 = MASKT_C0_3 >> 1;
  private static final long MASKT_C0_5 = MASKT_C0_4 >> 1;
  private static final long MASKT_C1 = MASK_C1;


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
    final long laa = (v0 & MASKT_A) << 18;
    final long lba =  v1 & MASKT_A;
    final long packedA = laa | lba;
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


//    System.err.println("TMASKAa=" + Utils.toBits(laa, NUMBER_BITS));
//    System.err.println("TMASKAb=" + Utils.toBits(lba, NUMBER_BITS));
//
//    System.err.println("TMASKBa=" + Utils.toBits(lab05, NUMBER_BITS));
//    System.err.println("TMASKBb=" + Utils.toBits(lab1, NUMBER_BITS));
//    System.err.println("TMASKBc=" + Utils.toBits(lbb05, NUMBER_BITS));
//    System.err.println("TMASKBd=" + Utils.toBits(lbb1, NUMBER_BITS));
//
//    System.err.println("TMASKCa=" + Utils.toBits(lac03, NUMBER_BITS));
//    System.err.println("TMASKCb=" + Utils.toBits(lac1, NUMBER_BITS));
//    System.err.println("TMASKCc=" + Utils.toBits(lbc03, NUMBER_BITS));
//    System.err.println("TMASKCd=" + Utils.toBits(lbc1, NUMBER_BITS));

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

