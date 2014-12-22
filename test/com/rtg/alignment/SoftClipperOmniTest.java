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
package com.rtg.alignment;

import junit.framework.TestCase;

/**
 */
public class SoftClipperOmniTest extends TestCase {

  public void testStartInsertClip() {
    final SoftClipperOmni clipper = new SoftClipperOmni(null, 5);
    int[] soft = clipper.softClipActions(ActionsHelper.build("==========", 0, 0), false);
    assertEquals(0, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("==========", ActionsHelper.toString(soft));

    soft = clipper.softClipActions(ActionsHelper.build("==I=======", 0, 2), false);
                                             //         -- -  new start is *2* because the insert doesn't exist on the template.
    assertEquals(2, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("SSS=======", ActionsHelper.toString(soft));

    soft = clipper.softClipActions(ActionsHelper.build("I=========", 0, 2), false);
    assertEquals(0, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("S=========", ActionsHelper.toString(soft));

    soft = clipper.softClipActions(ActionsHelper.build("=I=========", 0, 2), false);
    assertEquals(1, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("SS=========", ActionsHelper.toString(soft));

    soft = clipper.softClipActions(ActionsHelper.build("====I========", 0, 2), false);
    assertEquals(4, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("SSSSS========", ActionsHelper.toString(soft));

    soft = clipper.softClipActions(ActionsHelper.build("=====I======", 3, 2), false);
    assertEquals(3, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("=====I======", ActionsHelper.toString(soft));

    soft = clipper.softClipActions(ActionsHelper.build("====II======", 0, 3), false);
    assertEquals(4, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("SSSSSS======", ActionsHelper.toString(soft));

    soft = clipper.softClipActions(ActionsHelper.build("====III=====", 0, 3), false);
    assertEquals(4, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("SSSSSSS=====", ActionsHelper.toString(soft));
  }

  public void testStartDelClip() {
    final SoftClipperOmni clipper = new SoftClipperOmni(null, 5);
    int[] soft = clipper.softClipActions(ActionsHelper.build("==D=======", 0, 2), false);
                                                   //         ----  new start is *3*
    assertEquals(3, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("SS=======", ActionsHelper.toString(soft)); //but only 2 read nt are used up

    soft = clipper.softClipActions(ActionsHelper.build("D========", 0, 2), false); //ridiculous case that should never occur
    assertEquals(1, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("========", ActionsHelper.toString(soft));

    soft = clipper.softClipActions(ActionsHelper.build("====D========", 0, 2), false);
    assertEquals(5, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("SSSS========", ActionsHelper.toString(soft));

    soft = clipper.softClipActions(ActionsHelper.build("=====D======", 3, 2), false);
    assertEquals(3, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("=====D======", ActionsHelper.toString(soft));

    soft = clipper.softClipActions(ActionsHelper.build("====DD======", 0, 3), false);
    assertEquals(6, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("SSSS======", ActionsHelper.toString(soft));

    soft = clipper.softClipActions(ActionsHelper.build("====DDD=====", 0, 3), false);
    assertEquals(7, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("SSSS=====", ActionsHelper.toString(soft));
  }

  public void testEndInsertClip() {
    final SoftClipperOmni clipper = new SoftClipperOmni(null, 5);
    int[] soft = clipper.softClipActions(ActionsHelper.build("=======I==", 0, 2), false);
    assertEquals(0, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("=======SSS", ActionsHelper.toString(soft));

    soft = clipper.softClipActions(ActionsHelper.build("========I====", 0, 2), false);
    assertEquals(0, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("========SSSSS", ActionsHelper.toString(soft));

    soft = clipper.softClipActions(ActionsHelper.build("========I=====", 0, 2), false);
    assertEquals(0, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("========I=====", ActionsHelper.toString(soft));

    soft = clipper.softClipActions(ActionsHelper.build("======II====", 6, 2), false);
    assertEquals(6, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("======SSSSSS", ActionsHelper.toString(soft));

    soft = clipper.softClipActions(ActionsHelper.build("=====III====", 0, 3), false);
    assertEquals(0, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("=====SSSSSSS", ActionsHelper.toString(soft));
  }

  public void testEndInsertClipRc() {
    final SoftClipperOmni clipper = new SoftClipperOmni(null, 5);
    int[] soft = clipper.softClipActions(ActionsHelper.build("=======I==", 0, 2), true);
    assertEquals(2, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("=======SSS", ActionsHelper.toString(soft));

    soft = clipper.softClipActions(ActionsHelper.build("========I====", 0, 2), true);
    assertEquals(4, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("========SSSSS", ActionsHelper.toString(soft));

    soft = clipper.softClipActions(ActionsHelper.build("========I=====", 0, 2), true);
    assertEquals(0, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("========I=====", ActionsHelper.toString(soft));

    soft = clipper.softClipActions(ActionsHelper.build("======II====", 6, 2), true);
    assertEquals(10, ActionsHelper.zeroBasedTemplateStart(soft)); //6 + 4
    assertEquals("======SSSSSS", ActionsHelper.toString(soft));

    soft = clipper.softClipActions(ActionsHelper.build("=====III====", 0, 3), true);
    assertEquals(4, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("=====SSSSSSS", ActionsHelper.toString(soft));
  }

  public void testEndDelClip() {
    final SoftClipperOmni clipper = new SoftClipperOmni(null, 5);
    int[] soft = clipper.softClipActions(ActionsHelper.build("=======D==", 0, 2), false);
    assertEquals(0, ActionsHelper.zeroBasedTemplateStart(soft));

    soft = clipper.softClipActions(ActionsHelper.build("========D====", 0, 2), false);
    assertEquals(0, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("========SSSS", ActionsHelper.toString(soft));

    soft = clipper.softClipActions(ActionsHelper.build("========D=====", 0, 2), false);
    assertEquals(0, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("========D=====", ActionsHelper.toString(soft));

    soft = clipper.softClipActions(ActionsHelper.build("======DD====", 6, 2), false);
    assertEquals(6, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("======SSSS", ActionsHelper.toString(soft));

    soft = clipper.softClipActions(ActionsHelper.build("=====DDD====", 0, 3), false);
    assertEquals(0, ActionsHelper.zeroBasedTemplateStart(soft));
    assertEquals("=====SSSS", ActionsHelper.toString(soft));
  }

  public void testShortcuts() {
    final EditDistance ed = new SoftClipperOmni(new EditDistance() {
      @Override
      public int[] calculateEditDistance(byte[] read, int rlen, byte[] template, int zeroBasedStart, boolean rc, int maxScore, int maxShift, boolean cgLeft) {
        return ActionsHelper.build("=======D==", 4, 23);
      }
      @Override
      public void logStats() { }
    }, 0);

    final int[] actions = ed.calculateEditDistance(null, 10, null, 0, false, 90, 10, false);
    assertEquals(4, ActionsHelper.zeroBasedTemplateStart(actions));
    assertEquals("=======D==", ActionsHelper.toString(actions));

    final EditDistance ed2 = new SoftClipperOmni(new EditDistance() {
      @Override
      public int[] calculateEditDistance(byte[] read, int rlen, byte[] template, int zeroBasedStart, boolean rc, int maxScore, int maxShift, boolean cgLeft) {
        return ActionsHelper.build("=D=", 4, Integer.MAX_VALUE);
      }
      @Override
      public void logStats() { }
    }, 3);

    final int[] actions2 = ed2.calculateEditDistance(null, 10, null, 0, false, 90, 10, false);
    assertEquals(4, ActionsHelper.zeroBasedTemplateStart(actions2));
    assertEquals("=D=", ActionsHelper.toString(actions2));

  }

  public void testEndDel() {
    final EditDistance ed = new SoftClipperOmni(new EditDistance() {
      @Override
      public int[] calculateEditDistance(byte[] read, int rlen, byte[] template, int zeroBasedStart, boolean rc, int maxScore, int maxShift, boolean cgLeft) {
        return ActionsHelper.build("=================================================================================================DDD====", 0, 23);
      }
      @Override
      public void logStats() { }
    }, 5);

    final int[] actions = ed.calculateEditDistance(null, 101, null, 0, false, 90, 10, false);
    assertEquals(0, ActionsHelper.zeroBasedTemplateStart(actions));
    assertEquals("=================================================================================================SSSS", ActionsHelper.toString(actions));
  }

  public void testReverseComplimentDel() {
    final EditDistance ed = new SoftClipperOmni(new EditDistance() {
      @Override
      public int[] calculateEditDistance(byte[] read, int rlen, byte[] template, int zeroBasedStart, boolean rc, int maxScore, int maxShift, boolean cgLeft) {
        return ActionsHelper.build("=================================================================================================DDD====", 1, 23);
      }
      @Override
      public void logStats() { }
    }, 5);

    final int[] actions = ed.calculateEditDistance(null, 101, null, 1, true, 90, 10, false);

    assertEquals(8, ActionsHelper.zeroBasedTemplateStart(actions));
    assertEquals("=================================================================================================SSSS", ActionsHelper.toString(actions));
  }

  public void testBothEnds() {
    final SoftClipperOmni clipper = new SoftClipperOmni(null, 5);
    int[] soft = clipper.softClipActions(ActionsHelper.build("==I==============I==", 0, 0), false);
    assertEquals("SSS==============SSS", ActionsHelper.toString(soft));
    assertEquals(2, ActionsHelper.zeroBasedTemplateStart(soft));

    soft = clipper.softClipActions(ActionsHelper.build("==IIIIIIIIIIIIIIII==", 0, 0), false);
    assertEquals("==SSSSSSSSSSSSSSSSSS", ActionsHelper.toString(soft));
    assertEquals(0, ActionsHelper.zeroBasedTemplateStart(soft));
  }

  public void testClipExtension() {
    final SoftClipperOmni clipper = new SoftClipperOmni(null, 5);
    int[] soft = clipper.softClipActions(ActionsHelper.build("===IXXXXX=====", 0, 0), false);
    assertEquals("SSSSSSSSS=====", ActionsHelper.toString(soft));
    assertEquals(8, ActionsHelper.zeroBasedTemplateStart(soft));

    soft = clipper.softClipActions(ActionsHelper.build("===IXXXXXI====", 0, 0), false);
    assertEquals("===SSSSSSSSSSS", ActionsHelper.toString(soft));
    assertEquals(0, ActionsHelper.zeroBasedTemplateStart(soft));

    soft = clipper.softClipActions(ActionsHelper.build("===IXX==XXI===", 0, 0), false);
    assertEquals("SSSSSS==SSSSSS", ActionsHelper.toString(soft));
    assertEquals(5, ActionsHelper.zeroBasedTemplateStart(soft));

    soft = clipper.softClipActions(ActionsHelper.build("=IXX=====", 0, 0), false);
    assertEquals("SSSS=====", ActionsHelper.toString(soft));
    assertEquals(3, ActionsHelper.zeroBasedTemplateStart(soft));

    soft = clipper.softClipActions(ActionsHelper.build("=DXX=====", 0, 0), false);
    assertEquals("SSS=====", ActionsHelper.toString(soft));
    assertEquals(4, ActionsHelper.zeroBasedTemplateStart(soft));

    soft = clipper.softClipActions(ActionsHelper.build("=IXXXXXDXX=====", 0, 0), false);
    assertEquals("SSSSSSSSS=====", ActionsHelper.toString(soft));
    assertEquals(9, ActionsHelper.zeroBasedTemplateStart(soft));
  }
}
