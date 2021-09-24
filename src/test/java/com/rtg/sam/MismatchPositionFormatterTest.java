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
package com.rtg.sam;

import com.rtg.alignment.ActionsHelper;
import com.rtg.mode.DnaUtils;

import junit.framework.TestCase;

/**
 * Test class for {@link MismatchPositionFormatter}
 */
public class MismatchPositionFormatterTest extends TestCase {

  static final byte[] TEMPLATE = DnaUtils.encodeString("ACTGACTGACTGACTG");

  public void test() {
    assertEquals("7", MismatchPositionFormatter.actionsToMismatchPositions(ActionsHelper.build("=======", 0, 0), false, TEMPLATE));
    assertEquals("3G3", MismatchPositionFormatter.actionsToMismatchPositions(ActionsHelper.build("===X===", 0, 0), false, TEMPLATE));
    assertEquals("0A5T0", MismatchPositionFormatter.actionsToMismatchPositions(ActionsHelper.build("X=====X", 0, 0), false, TEMPLATE));
    assertEquals("0A0^C0T5", MismatchPositionFormatter.actionsToMismatchPositions(ActionsHelper.build("XDX=====", 0, 0), false, TEMPLATE));
    assertEquals("0^A5^T0", MismatchPositionFormatter.actionsToMismatchPositions(ActionsHelper.build("D=====D", 0, 0), false, TEMPLATE));
    assertEquals("4A2", MismatchPositionFormatter.actionsToMismatchPositions(ActionsHelper.build("==X====", 0, 0), true, TEMPLATE));
    assertEquals("2T0G3", MismatchPositionFormatter.actionsToMismatchPositions(ActionsHelper.build("===XX==", 0, 0), true, TEMPLATE));
    assertEquals("1^C5", MismatchPositionFormatter.actionsToMismatchPositions(ActionsHelper.build("=D=====", 0, 0), false, TEMPLATE));
    assertEquals("2T1^AC1", MismatchPositionFormatter.actionsToMismatchPositions(ActionsHelper.build("=DD=X==", 0, 0), true, TEMPLATE));
    assertEquals("2G1^CT1", MismatchPositionFormatter.actionsToMismatchPositions(ActionsHelper.build("=DD=X==", 1, 0), true, TEMPLATE));
    assertEquals("1C1^GA1", MismatchPositionFormatter.actionsToMismatchPositions(ActionsHelper.build("=DD=X==", -1, 0), true, TEMPLATE));
    assertEquals("5", MismatchPositionFormatter.actionsToMismatchPositions(ActionsHelper.build("======", 11, 0), false, TEMPLATE));
    assertEquals("4C0", MismatchPositionFormatter.actionsToMismatchPositions(ActionsHelper.build("==N=II=BX", 0, 0), false, TEMPLATE));
  }

}
