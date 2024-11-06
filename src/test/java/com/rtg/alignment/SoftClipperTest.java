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
package com.rtg.alignment;

import junit.framework.TestCase;

/**
 */
public class SoftClipperTest extends TestCase {

  public void testStartInsertClip() {
    final SoftClipper clipper = new SoftClipper(null, 5, 0);
    int[] soft = clipper.softClipActions(true, true, ActionsHelper.build("==========", 0, 0));
    assertEquals(0, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("==========", ActionsHelper.toString(soft));

    soft = clipper.softClipActions(true, true, ActionsHelper.build("==I=======", 0, 2));
                                             //         -- -  new start is *2* because the insert doesn't exist on the template.
    assertEquals(2, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("SSS=======", ActionsHelper.toString(soft));

    soft = clipper.softClipActions(true, true, ActionsHelper.build("I=========", 0, 2));
    assertEquals(0, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("S=========", ActionsHelper.toString(soft));

    soft = clipper.softClipActions(true, true, ActionsHelper.build("=I=========", 0, 2));
    assertEquals(1, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("SS=========", ActionsHelper.toString(soft));

    soft = clipper.softClipActions(true, true, ActionsHelper.build("====I========", 0, 2));
    assertEquals(4, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("SSSSS========", ActionsHelper.toString(soft));

    soft = clipper.softClipActions(true, true, ActionsHelper.build("=====I======", 3, 2));
    assertEquals(3, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("=====I======", ActionsHelper.toString(soft));

    soft = clipper.softClipActions(true, true, ActionsHelper.build("====II======", 0, 3));
    assertEquals(4, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("SSSSSS======", ActionsHelper.toString(soft));

    soft = clipper.softClipActions(true, true, ActionsHelper.build("====III=====", 0, 3));
    assertEquals(4, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("SSSSSSS=====", ActionsHelper.toString(soft));
  }

  public void testStartDelClip() {
    final SoftClipper clipper = new SoftClipper(null, 5, 0);
    int[] soft = clipper.softClipActions(true, true, ActionsHelper.build("==D=======", 0, 2));
                                                   //         ----  new start is *3*
    assertEquals(3, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("SS=======", ActionsHelper.toString(soft)); //but only 2 read nt are used up

    soft = clipper.softClipActions(true, true, ActionsHelper.build("D========", 0, 2)); //ridiculous case that should never occur
    assertEquals(1, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("========", ActionsHelper.toString(soft));

    soft = clipper.softClipActions(true, true, ActionsHelper.build("====D========", 0, 2));
    assertEquals(5, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("SSSS========", ActionsHelper.toString(soft));

    soft = clipper.softClipActions(true, true, ActionsHelper.build("=====D======", 3, 2));
    assertEquals(3, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("=====D======", ActionsHelper.toString(soft));

    soft = clipper.softClipActions(true, true, ActionsHelper.build("====DD======", 0, 3));
    assertEquals(6, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("SSSS======", ActionsHelper.toString(soft));

    soft = clipper.softClipActions(true, true, ActionsHelper.build("====DDD=====", 0, 3));
    assertEquals(7, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("SSSS=====", ActionsHelper.toString(soft));
  }

  public void testEndInsertClip() {
    final SoftClipper clipper = new SoftClipper(null, 5, 0);
    int[] soft = clipper.softClipActions(true, true, ActionsHelper.build("=======I==", 0, 2));
    assertEquals(0, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("=======SSS", ActionsHelper.toString(soft));

    soft = clipper.softClipActions(true, true, ActionsHelper.build("========I====", 0, 2));
    assertEquals(0, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("========SSSSS", ActionsHelper.toString(soft));

    soft = clipper.softClipActions(true, true, ActionsHelper.build("========I=====", 0, 2));
    assertEquals(0, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("========I=====", ActionsHelper.toString(soft));

    soft = clipper.softClipActions(true, true, ActionsHelper.build("======II====", 6, 2));
    assertEquals(6, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("======SSSSSS", ActionsHelper.toString(soft));

    soft = clipper.softClipActions(true, true, ActionsHelper.build("=====III====", 0, 3));
    assertEquals(0, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("=====SSSSSSS", ActionsHelper.toString(soft));
  }

  public void testEndDelClip() {
    final SoftClipper clipper = new SoftClipper(null, 5, 0);
    int[] soft = clipper.softClipActions(true, true, ActionsHelper.build("=======D==", 0, 2));
    assertEquals(0, ActionsHelper.zeroBasedTemplateStart(soft));

    soft = clipper.softClipActions(true, true, ActionsHelper.build("========D====", 0, 2));
    assertEquals(0, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("========SSSS", ActionsHelper.toString(soft));

    soft = clipper.softClipActions(true, true, ActionsHelper.build("========D=====", 0, 2));
    assertEquals(0, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("========D=====", ActionsHelper.toString(soft));

    soft = clipper.softClipActions(true, true, ActionsHelper.build("======DD====", 6, 2));
    assertEquals(6, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("======SSSS", ActionsHelper.toString(soft));

    soft = clipper.softClipActions(true, true, ActionsHelper.build("=====DDD====", 0, 3));
    assertEquals(0, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("=====SSSS", ActionsHelper.toString(soft));
  }

  private abstract class UedAdaptor implements UnidirectionalEditDistance {
    @Override
    public void logStats() { }
    @Override
    public int[] calculateEditDistanceFixedStart(byte[] read, int readStartPos, int readEndPos, byte[] template, int templateStartPos, int maxScore, int maxShift) {
      return null;
    }
    @Override
    public int[] calculateEditDistanceFixedEnd(byte[] read, int readStartPos, int readEndPos, byte[] template, int templateExpectedStartPos, int templateEndPos, int maxScore, int maxShift) {
      return null;
    }
    @Override
    public int[] calculateEditDistanceFixedBoth(byte[] read, int readStartPos, int readEndPos, byte[] template, int templateStartPos, int templateEndPos, int maxScore, int maxShift) {
      return null;
    }
  }
  public void testShortcuts() {
    final UnidirectionalEditDistance ed = new SoftClipper(new UedAdaptor() {
      @Override
      public int[] calculateEditDistance(byte[] read, int rlen, byte[] template, int zeroBasedStart, int maxScore, int maxShift, boolean cgLeft) {
        return ActionsHelper.build("=======D==", 4, 23);
      }

    }, 0, 0);

    final int[] actions = ed.calculateEditDistance(null, 10, null, 0, 90, 10, false);
    assertEquals(4, ActionsHelper.zeroBasedTemplateStart(actions));
    assertEquals("=======D==", ActionsHelper.toString(actions));

    final UnidirectionalEditDistance ed2 = new SoftClipper(new UedAdaptor() {
      @Override
      public int[] calculateEditDistance(byte[] read, int rlen, byte[] template, int zeroBasedStart, int maxScore, int maxShift, boolean cgLeft) {
        return ActionsHelper.build("=D=", 4, Integer.MAX_VALUE);
      }
    }, 3, 0);

    final int[] actions2 = ed2.calculateEditDistance(null, 10, null, 0, 90, 10, false);
    assertEquals(4, ActionsHelper.zeroBasedTemplateStart(actions2));
    assertEquals("=D=", ActionsHelper.toString(actions2));

  }

  public void testEndDel() {
    final UnidirectionalEditDistance ed = new SoftClipper(new UedAdaptor() {
      @Override
      public int[] calculateEditDistance(byte[] read, int rlen, byte[] template, int zeroBasedStart, int maxScore, int maxShift, boolean cgLeft) {
        return ActionsHelper.build("=================================================================================================DDD====", 0, 23);
      }
    }, 5, 0);

    final int[] actions = ed.calculateEditDistance(null, 101, null, 0, 90, 10, false);
    assertEquals(0, ActionsHelper.zeroBasedTemplateStart(actions));
    assertEquals("=================================================================================================SSSS", ActionsHelper.toString(actions));
  }

  public void testBothEnds() {
    final SoftClipper clipper = new SoftClipper(null, 5, 0);
    int[] soft = clipper.softClipActions(true, true, ActionsHelper.build("==I==============I==", 0, 0));
    assertEquals("SSS==============SSS", ActionsHelper.toString(soft));
    assertEquals(2, ActionsHelper.zeroBasedTemplateStart(soft));

    soft = clipper.softClipActions(true, true, ActionsHelper.build("==IIIIIIIIIIIIIIII==", 0, 0));
    assertEquals("==SSSSSSSSSSSSSSSSSS", ActionsHelper.toString(soft));
    assertEquals(0, ActionsHelper.zeroBasedTemplateStart(soft));
  }

  public void testClipExtension() {
    final SoftClipper clipper = new SoftClipper(null, 5, 0);
    int[] soft = clipper.softClipActions(true, true, ActionsHelper.build("===IXXXXX=====", 0, 0));
    assertEquals("SSSSSSSSS=====", ActionsHelper.toString(soft));
    assertEquals(8, ActionsHelper.zeroBasedTemplateStart(soft));

    soft = clipper.softClipActions(true, true, ActionsHelper.build("===IXXXXXI====", 0, 0));
    assertEquals("===SSSSSSSSSSS", ActionsHelper.toString(soft));
    assertEquals(0, ActionsHelper.zeroBasedTemplateStart(soft));

    soft = clipper.softClipActions(true, true, ActionsHelper.build("===IXX==XXI===", 0, 0));
    assertEquals("SSSSSS==SSSSSS", ActionsHelper.toString(soft));
    assertEquals(5, ActionsHelper.zeroBasedTemplateStart(soft));

    soft = clipper.softClipActions(true, true, ActionsHelper.build("=IXX=====", 0, 0));
    assertEquals("SSSS=====", ActionsHelper.toString(soft));
    assertEquals(3, ActionsHelper.zeroBasedTemplateStart(soft));

    soft = clipper.softClipActions(true, true, ActionsHelper.build("=DXX=====", 0, 0));
    assertEquals("SSS=====", ActionsHelper.toString(soft));
    assertEquals(4, ActionsHelper.zeroBasedTemplateStart(soft));

    soft = clipper.softClipActions(true, true, ActionsHelper.build("=IXXXXXDXX=====", 0, 0));
    assertEquals("SSSSSSSSS=====", ActionsHelper.toString(soft));
    assertEquals(9, ActionsHelper.zeroBasedTemplateStart(soft));
  }

  public void testClipMismatch() {
    final SoftClipper clipper = new SoftClipper(null, 0, 3);
    int[] soft = clipper.softClipActions(true, true, ActionsHelper.build("==============", 0, 0));
    assertEquals("==============", ActionsHelper.toString(soft));
    assertEquals(0, ActionsHelper.zeroBasedTemplateStart(soft));

    // First mismatch before enough anchoring bases seen
    soft = clipper.softClipActions(true, true, ActionsHelper.build("=XXXXXXXX=====", 0, 0));
    assertEquals("SSSSSSSSS=====", ActionsHelper.toString(soft));
    assertEquals(9, ActionsHelper.zeroBasedTemplateStart(soft));

    // First mismatch before enough anchoring bases seen
    soft = clipper.softClipActions(true, true, ActionsHelper.build("==XXXXXXX=====", 0, 0));
    assertEquals("SSSSSSSSS=====", ActionsHelper.toString(soft));
    assertEquals(9, ActionsHelper.zeroBasedTemplateStart(soft));

    // First mismatch only after enough anchoring bases seen
    soft = clipper.softClipActions(true, true, ActionsHelper.build("===XXXXXX=====", 0, 0));
    assertEquals("===XXXXXX=====", ActionsHelper.toString(soft));
    assertEquals(0, ActionsHelper.zeroBasedTemplateStart(soft));

    soft = clipper.softClipActions(true, true, ActionsHelper.build("==X=X=========", 0, 0));
    assertEquals("SSS=X=========", ActionsHelper.toString(soft));
    assertEquals(3, ActionsHelper.zeroBasedTemplateStart(soft));

    soft = clipper.softClipActions(true, true, ActionsHelper.build("==XIIII=======", 0, 0));
    assertEquals("SSSSSSS=======", ActionsHelper.toString(soft));
    assertEquals(3, ActionsHelper.zeroBasedTemplateStart(soft));

    soft = clipper.softClipActions(true, true, ActionsHelper.build("==XDDDD=======", 0, 0));
    assertEquals("SSS=======", ActionsHelper.toString(soft));
    assertEquals(7, ActionsHelper.zeroBasedTemplateStart(soft));

    soft = clipper.softClipActions(true, true, ActionsHelper.build("X==", 0, 0));
    assertEquals("SSS", ActionsHelper.toString(soft));
    assertEquals(0, ActionsHelper.zeroBasedTemplateStart(soft));

    soft = clipper.softClipActions(true, true, ActionsHelper.build("==X", 0, 0));
    assertEquals("==S", ActionsHelper.toString(soft));
    assertEquals(0, ActionsHelper.zeroBasedTemplateStart(soft));
  }

  public void testClipMinMatch() {
    final SoftClipper clipper = new SoftClipper(null, 5, 5, 5);
    final int[] actions = ActionsHelper.build("X=IIIIIIIIIIIIIIII=X", 0, 0);
    assertEquals(0, ActionsHelper.alignmentScore(actions));
    final int[] soft = clipper.softClipActions(true, true, actions);
    assertEquals("S=SSSSSSSSSSSSSSSSSS", ActionsHelper.toString(soft));
    assertEquals(1, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals(Integer.MAX_VALUE, ActionsHelper.alignmentScore(soft));
  }
}
