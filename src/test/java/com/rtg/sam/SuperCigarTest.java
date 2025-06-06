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

package com.rtg.sam;

import com.rtg.alignment.ActionsHelper;
import com.rtg.mode.DNA;
import com.rtg.mode.DnaUtils;

import junit.framework.TestCase;

/**
 */
public class SuperCigarTest extends TestCase {

  public void testActionsToSuperCigarFwd() {
    assertEquals("", SuperCigar.actionsToSuperCigar(ActionsHelper.build("", 0, 0), false, 0));
    assertEquals("1X", SuperCigar.actionsToSuperCigar(ActionsHelper.build("X", 0, 0), false, 1));
    assertEquals("1=2X2=3B1=3I1X2=4N3=2D1=", SuperCigar.actionsToSuperCigar(ActionsHelper.build("=XX==BBB=IIIX==NNNN===DD=", 0, 0), false, 50));
  }

  public void testActionsToSuperCigarRev() {
    assertEquals("", SuperCigar.actionsToSuperCigar(ActionsHelper.build("", 0, 0), true, 50));
    assertEquals("7=", SuperCigar.actionsToSuperCigar(ActionsHelper.build("=======", 0, 0), true, 50));
    assertEquals("1=2D3=4N2=1X3I1=3B2=2X1=", SuperCigar.actionsToSuperCigar(ActionsHelper.build("=XX==BBB=IIIX==NNNN===DD=", 0, 0), true, 50));
  }

  public void testReadDelta() {
    assertEquals("", SuperCigar.readDelta(ActionsHelper.build("", 0, 0), new byte[] {}, 0));
    assertEquals("NACGT", SuperCigar.readDelta(ActionsHelper.build("XXXXX", 0, 0), new byte[] {0, 1, 2, 3, 4}, 5));
    final String read = "CAGAAGCTCGAACCCG";
    assertEquals("CGCTCG", SuperCigar.readDelta(ActionsHelper.build("X=X==BBB=IIIX==NNNN===DD=", 0, 0), DnaUtils.encodeString(read), 50));
  }

  public void testSoftClip() {
    assertEquals("2S", SuperCigar.actionsToSuperCigar(ActionsHelper.build("==", -2, 0), false, 50));
    assertEquals("1S1=", SuperCigar.actionsToSuperCigar(ActionsHelper.build("==", -1, 0), false, 50));
    assertEquals("2S2=", SuperCigar.actionsToSuperCigar(ActionsHelper.build("====", -2, 0), false, 50));
    assertEquals("5S2B2S2=", SuperCigar.actionsToSuperCigar(ActionsHelper.build("=====BB====", -5, 0), false, 50));
    assertEquals("5S2B4S1=", SuperCigar.actionsToSuperCigar(ActionsHelper.build("=====BB=====", -7, 0), false, 50));
  }
  public void testSoftDelta() {
    assertEquals("NA", SuperCigar.readDelta(ActionsHelper.build("==", -2, 0), new byte[] {0, 1}, 50));
    assertEquals("A", SuperCigar.readDelta(ActionsHelper.build("==", -1, 0), new byte[] {1, 2}, 50));
    assertEquals("AC", SuperCigar.readDelta(ActionsHelper.build("====", -2, 0), new byte[] {1, 2, 3, 4}, 50));
    assertEquals("ACGTACG", SuperCigar.readDelta(ActionsHelper.build("=====BB====", -5, 0), DnaUtils.encodeString("ACGTACGTA"), 50));
    assertEquals("ACGTACGTA", SuperCigar.readDelta(ActionsHelper.build("=====BB=====", -7, 0), DnaUtils.encodeString("ACGTACGTAC"), 50));
  }

  public void testRealWorldSoftClip() {
   int[] actions = ActionsHelper.build("=====BB====================NNNNN=========X", -3, 2);
   assertEquals("3S2=2B20=5N9=1X", SuperCigar.actionsToSuperCigar(actions, false, 50));
   assertEquals("AGCA", SuperCigar.readDelta(actions, DnaUtils.encodeString("AGCCCCCACACGTAAATAAGACATC     ACGATGATCA".replaceAll(" ", "")), 36));
   assertEquals("5S2B2S8=2N10=7N10=", SuperCigar.actionsToSuperCigar(ActionsHelper.build("=====BB==========NN==========NNNNNNN==========", -5, 3), false, 50));

   assertEquals("10=7N10=2N1X6=3S2B5S", SuperCigar.actionsToSuperCigar(ActionsHelper.build("==========NNNNNNN==========NNX=========BB=====", 0, 5), false, 36));
   assertEquals("10=6N18=2S2B5S", SuperCigar.actionsToSuperCigar(ActionsHelper.build("==========NNNNNN====================BB=====", 0, 5), false, 34));

   actions = ActionsHelper.build("==========NNNNNNN====================BB=====", 0, 5);
   assertEquals("10=7N19=1S2B1=4S", SuperCigar.actionsToSuperCigar(actions, false, 36));
   assertEquals("AATCA", SuperCigar.readDelta(actions, DnaUtils.encodeString("AGCCCACACG     TAAATAAGACATCACGATGAGATCA".replaceAll(" ", "")), 36));
  }

  public void testNs() {
    final int[] actions = ActionsHelper.build("===XX===", 0, 0);
    assertEquals("3=2X3=", SuperCigar.actionsToSuperCigar(actions, false, 50));
    assertEquals("NA", SuperCigar.readDelta(actions, DnaUtils.encodeString("GGGNAGGG"), 50));
  }

  public void testUnknownDeltas() {
    final int[] actions = ActionsHelper.build("===TR===", 0, 0);
    assertEquals("3=1T1R3=", SuperCigar.actionsToSuperCigar(actions, false, 50));
    assertEquals("T", SuperCigar.readDelta(actions, DnaUtils.encodeString("ACGTGCAT"), 50));
  }

  public void testReverseSoftClip() {
    int[] actions = ActionsHelper.build("X====BB====================NNNNNN==========", 0, 1);
    assertEquals("10=6N20=2B4=1S", SuperCigar.actionsToSuperCigar(actions, true, "TCTTTCTCAGGGATGTTCTTTGCTGAGAAAAAGAATTC".length()));

    actions = ActionsHelper.build("==========NNNNNNN====================BB==XXX", -3, 3);
    assertEquals("3S2=2B20=7N10=", SuperCigar.actionsToSuperCigar(actions, true, 100));
  }

  public void testEndTemplateProblems() {
    final int[] actions = ActionsHelper.build("XXXXXBBXXXXXXXX============NNNNNN=======X==", 728, 3);
    assertEquals("ACGTACGTACGTAACGTACGTACG", SuperCigar.readDelta(actions, DNA.stringDNAtoByte("ACGTACGTACGTACGTACGTACGTACGTACGTACG"), 750));
  }
}
