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
package com.rtg.util;

import com.rtg.mode.BidirectionalFrame;

import junit.framework.TestCase;

/**
 *
 *
 */
public class EnumHelperTest extends TestCase {

  /**
   */
  public EnumHelperTest(final String name) {
    super(name);
  }

  public void testEnum() {
    EnumHelper<BidirectionalFrame> helper = new EnumHelper<>(BidirectionalFrame.class, new BidirectionalFrame[] {BidirectionalFrame.FORWARD, BidirectionalFrame.REVERSE});
    try {
      helper.valueOf("temp");
      fail();
    } catch (final IllegalArgumentException iae) {
      assertEquals("temp is not a valid enum value of type: " + BidirectionalFrame.class.toString(), iae.getMessage());
    }

    assertEquals(BidirectionalFrame.FORWARD, helper.valueOf("FORWARD"));
    assertEquals(BidirectionalFrame.REVERSE, helper.valueOf("REVERSE"));

    assertEquals(2, helper.values().length);
    assertEquals(BidirectionalFrame.FORWARD, helper.values()[0]);
    assertEquals(BidirectionalFrame.REVERSE, helper.values()[1]);

    assertEquals(2, helper.names().length);
    assertEquals("FORWARD", helper.names()[0]);
    assertEquals("REVERSE", helper.names()[1]);
  }



}
