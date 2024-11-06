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

package com.rtg.visualization;

import java.util.HashMap;

import com.rtg.mode.DnaUtils;
import com.rtg.sam.BadSuperCigarException;

import htsjdk.samtools.SAMRecord;

import junit.framework.TestCase;

/**
 */
public class SuperCigarUnrollerTest extends TestCase {



  public void testSuperCigarUnrolling() throws BadSuperCigarException {
    // 694 179 simulatedSequence1 19 255 21=1X1=6N10= = 224 205
    // AACGGATTTAGTGCCCCTGCTTAAGTTTACACA $\EC3(NP;M'JS>IB;IXL[\!MDRR_.SOD=
    // AS:i:1 MQ:i:255
    // 0123456789012345678901234567890123456789
    // AACGGATTTAGTGCCCCTGCTCANNNNNNAGTTTACACA
    // AACGG
    // GGATTTAGTGCCCCTGCTTA AGTTTACACA
    // XU:Z:5=2B18=1X1=6N10= XQ:Z:R: XR:Z:T XA:i:4 IH:i:1

    final SAMRecord rec = SamTestHelper.getSuperCigarSAMRecord("21=1X1=6N10=", "AACGGATTTAGTGCCCCTGCTTAAGTTTACACA", "5=2B18=1X1=6N10=", "R:", "T");

    final SuperCigarUnroller scu = new SuperCigarUnroller();
    final HashMap<Integer, String> map = scu.unroll(rec, DnaUtils.encodeString("AACGGATTTAGTGCCCCTGCTCANNNNNNAGTTTACACA"));
    //printMap(map);
    assertEquals(2, map.size());
    assertEquals("AACGG", map.get(0));
    assertEquals("GGATTTAGTGCCCCTGCTTA      AGTTTACACA", map.get(3));
  }


  public void testLargerOverlap() throws BadSuperCigarException {
 // 1070 131 simulatedSequence1 1256 255 10=6N7=2X13= = 1059 -197
    //
    // GGTCGCGACG      GAAGTCCGGATTCAGTATAGTC L2,D.[F?426!U-290Z\%K6T^-IF*"@D: AS:i:2
    // GGTCGCGACGNNNNNNGAAGTCCTTATTCAGTATAG
    //                                  TAGTC
    // MQ:i:255

    // XU:Z:10=6N7=2X11=3B5= XQ:Z:^*- XR:Z:GG XA:i:4 IH:i:1
    final SAMRecord rec2 = SamTestHelper.getSuperCigarSAMRecord("10=6N7=2X13=", "GGTCGCGACGGAAGTCCGGATTCAGTATAGTC", "10=6N7=2X11=3B5=", "^*-", "GG");

    final SuperCigarUnroller scu2 = new SuperCigarUnroller();
    final HashMap<Integer, String> map2 = scu2.unroll(rec2, DnaUtils.encodeString("GGTCGCGACGNNNNNNGAAGTCCTTATTCAGTATAGTC"));
    //printMap(map2);
    assertEquals(2, map2.size());
  }

//  private void printMap(HashMap<Integer, String> map) {
//    for (Map.Entry<Integer, String> e : map.entrySet()) {
//      System.err.println(getSpace(e.getKey()) + " " + e.getValue());
//    }
//  }
//
//  private String getSpace(int key) {
//    final StringBuilder sb = new StringBuilder();
//    for (int i = 0; i < key; ++i) {
//      sb.append(' ');
//    }
//    return sb.toString();
//  }


}
